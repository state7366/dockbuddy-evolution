# -*- coding: utf-8 -*-
"""libdock.ranking —— 多目标突变候选排序层(B类, #311)。

================================================================================
定位与诚实声明(最高优先级, 见 verdict):
- 本层是**实验优先级启发式**, 不是活性/功能预测器。它把若干"已算好、各自带自身
  verdict 免责"的指标(ΔΔG、取向占比、空腔变化、上位性、共识禁忌命中…)统一做多
  目标排序, 仅用于"先算/先验证哪几个候选"的**排程建议**, 绝不产出"最优突变 =
  活性最高"的结论。
- 归一化口径: 逐目标 min-max 归一化到 [0,1], 按 objective.direction 翻转使
  "1 = 更优"; 加权求和得综合分(weight 缺省 1.0)。缺失值默认中性 0.5(不帮不害),
  可由 on_missing 改为 "worst"(0.0) 或 "exclude"(整项退出排序)。
- 帕累托口径: 在所有目标"越大越优"的归一化空间里求非支配前沿; 前沿 pareto_rank=0。
- 硬门禁(exclude_if / blocker_keys): 任一硬门失败(结构 QC 不过、撞共识禁忌…)的
  候选被**剔除出排序**(单列 excluded), 不参与加权/前沿, 且不得与通过者混排。
- verdict 强制禁写: "排序第一 = 活性最高/功能最优"、"加权分高 = 实验必成功"、
  "帕累托前沿 = 一定可行"; 可写: "排序仅为实验优先级建议"、"结论须湿实验收口"。
================================================================================
"""
from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple


def _resolve_metric(item: Dict[str, Any], key: str) -> Any:
    """指标可在 item 顶层或 item["metrics"] 子 dict 里。"""
    if key in item:
        return item[key]
    m = item.get("metrics")
    if isinstance(m, dict) and key in m:
        return m[key]
    return None


def make_blocker_predicate(require_true: Optional[List[str]] = None,
                           require_false: Optional[List[str]] = None
                           ) -> Callable[[Dict[str, Any]], Optional[str]]:
    """构造 exclude_predicate:
      - require_true 中任一指标为 falsy(0/False/None/"")→ 排除;
      - require_false 中任一指标为 truthy → 排除。
    返回谓词: 排除时给原因字符串, 否则 None。"""
    rt = list(require_true or [])
    rf = list(require_false or [])

    def _pred(item: Dict[str, Any]) -> Optional[str]:
        for k in rt:
            v = _resolve_metric(item, k)
            if not v:
                return f"blocker require_true 未过: {k}={v!r}"
        for k in rf:
            v = _resolve_metric(item, k)
            if v:
                return f"blocker require_false 触发: {k}={v!r}"
        return None

    return _pred


def ranking_layer(items: List[Dict[str, Any]], *,
                  objectives: List[Dict[str, Any]],
                  weights: Optional[Dict[str, float]] = None,
                  method: str = "both",
                  exclude_predicate: Optional[Callable[[Dict[str, Any]], Optional[str]]] = None,
                  on_missing: str = "neutral",
                  blocker_keys: Optional[List[str]] = None
                  ) -> Dict[str, Any]:
    """多目标排序主入口。

    参数:
      items       候选列表, 每项至少含 "id"(唯一) 与各 objective.key 对应指标。
                  指标可放顶层或 items[i]["metrics"] 子 dict。
      objectives  目标定义列表, 每项 {name, key, direction("min"/"max", 默认 "max"),
                  weight(默认 1.0)}。direction="min" 表示原值越小越优(如 ΔΔG)。
      weights     可选 {name: w} 覆盖 objectives 里的 weight。
      method      "weighted" | "pareto" | "both"(默认, 输出含两者)。
      exclude_predicate  自定义谓词(item)->Optional[str]; 返回原因即剔除。
      on_missing  "neutral"(默认, 缺失目标记 0.5) | "worst"(0.0) | "exclude"(整项退出)。
      blocker_keys  便捷硬门: 这些指标为 falsy 即剔除(等价于 make_blocker_predicate
                  (require_true=blocker_keys))。

    返回: {ok, method, n_active, n_excluded, ranked[], excluded[], objectives_stats,
           verdict}
    ranked[i] = {rank, id, score, pareto_rank, components{name:norm}, incomplete, item}
    """
    if method not in ("weighted", "pareto", "both"):
        return {"ok": False, "error": f"未知 method: {method!r}",
                "ranked": [], "excluded": [], "verdict": _verdict([], reason="未知 method")}
    if not isinstance(items, list) or not items:
        return {"ok": False, "error": "items 为空或非 list", "ranked": [],
                "excluded": [], "verdict": _verdict([], reason="空输入")}
    if not objectives:
        return {"ok": False, "error": "objectives 为空", "ranked": [],
                "excluded": [], "verdict": _verdict([], reason="缺 objectives")}

    # 归一化目标定义(权重优先取 weights 覆盖)
    objs: List[Dict[str, Any]] = []
    for o in objectives:
        nm = o["name"]
        w = float((weights or {}).get(nm, o.get("weight", 1.0)))
        objs.append({"name": nm, "key": o["key"],
                     "direction": o.get("direction", "max"), "weight": w})

    # 合并排除谓词
    preds: List[Callable[[Dict[str, Any]], Optional[str]]] = []
    if blocker_keys:
        preds.append(make_blocker_predicate(require_true=list(blocker_keys)))
    if exclude_predicate:
        preds.append(exclude_predicate)

    def _excl(it: Dict[str, Any]) -> Optional[str]:
        for p in preds:
            r = p(it)
            if r:
                return r
        return None

    # 1) 分离硬门失败的候选
    excluded: List[Dict[str, Any]] = []
    active: List[Dict[str, Any]] = []
    for it in items:
        iid = it.get("id", f"item{len(active)}")
        why = _excl(it)
        if why:
            excluded.append({"id": iid, "reason": why, "item": it})
        else:
            active.append(it)

    # 2) 逐目标 min-max 归一化(只用 active 中有值的候选)
    norm: Dict[str, Dict[Any, float]] = {}
    stats: Dict[str, Any] = {}
    for o in objs:
        key, direction = o["key"], o["direction"]
        raw: Dict[Any, float] = {}
        for it in active:
            v = _resolve_metric(it, key)
            if isinstance(v, (int, float)) and not isinstance(v, bool):
                raw[it.get("id")] = float(v)
        if not raw:
            stats[o["name"]] = {"min": None, "max": None, "present": 0, "direction": direction}
            norm[o["name"]] = {}
            continue
        lo, hi = min(raw.values()), max(raw.values())
        rng = (hi - lo) or 1.0
        per: Dict[Any, float] = {}
        for iid, val in raw.items():
            n = (val - lo) / rng            # 0..1, 越大=原值越大
            if direction == "min":
                n = 1.0 - n                 # 翻转: 越小越优 → 1
            per[iid] = n
        stats[o["name"]] = {"min": lo, "max": hi, "present": len(raw),
                            "direction": direction}
        norm[o["name"]] = per

    wsum = sum(o["weight"] for o in objs) or 1.0

    # 3) 计算加权分(缺失值按 on_missing 处理)
    ranked: List[Dict[str, Any]] = []
    for it in active:
        iid = it.get("id")
        comp: Dict[str, Optional[float]] = {}
        raw_metrics: Dict[str, Any] = {}
        for o in objs:
            nval = norm[o["name"]].get(iid)
            comp[o["name"]] = (None if nval is None else round(nval, 4))
            raw_metrics[o["key"]] = _resolve_metric(it, o["key"])
        if on_missing == "exclude" and any(v is None for v in comp.values()):
            missing = [nm for nm, v in comp.items() if v is None]
            excluded.append({"id": iid, "reason": f"目标缺失: {missing}", "item": it})
            continue
        score = 0.0
        for o in objs:
            nval = norm[o["name"]].get(iid)
            if nval is None:
                nval = 0.0 if on_missing == "worst" else 0.5
            score += o["weight"] * nval
        score = score / wsum
        ranked.append({
            "id": iid, "score": round(score, 4),
            "components": {k: (round(v, 4) if v is not None else None)
                          for k, v in comp.items()},
            "metrics": raw_metrics,
            "incomplete": any(v is None for v in comp.values()),
            "_norm": {nm: (norm[nm].get(iid) if norm[nm].get(iid) is not None
                           else (0.0 if on_missing == "worst" else 0.5))
                      for nm in norm},
            "item": it,
        })

    # 4) 帕累托前沿(在所有目标"越大越优"的归一化空间求非支配集)
    pareto_rank: Dict[Any, int] = {}
    if method in ("pareto", "both"):
        pts: List[Tuple[Any, List[float]]] = [
            (r["id"], [r["_norm"][o["name"]] for o in objs]) for r in ranked]
        remaining = list(pts)
        pr = 0
        while remaining:
            nondom = []
            for i, (idi, veci) in enumerate(remaining):
                dominated = False
                for j, (idj, vecj) in enumerate(remaining):
                    if i == j:
                        continue
                    if _dominates(vecj, veci):
                        dominated = True
                        break
                if not dominated:
                    nondom.append((idi, veci))
            for idi, _ in nondom:
                pareto_rank[idi] = pr
            remaining = [(i, v) for (i, v) in remaining if pareto_rank.get(i, -1) < 0]
            pr += 1

    # 5) 排序 + 附 rank
    if method == "weighted":
        ranked.sort(key=lambda r: -r["score"])
    elif method == "pareto":
        ranked.sort(key=lambda r: (pareto_rank.get(r["id"], 10 ** 9), -r["score"]))
    else:  # both: 前沿优先, 同前沿按加权分降序
        ranked.sort(key=lambda r: (pareto_rank.get(r["id"], 10 ** 9), -r["score"]))

    for idx, r in enumerate(ranked, 1):
        r["rank"] = idx
        r["pareto_rank"] = pareto_rank.get(r["id"], None)
        r.pop("_norm", None)

    top_id = ranked[0]["id"] if ranked else None
    return {
        "ok": True, "method": method,
        "n_active": len(ranked), "n_excluded": len(excluded),
        "n_objectives": len(objs),
        "objectives": [{"name": o["name"], "key": o["key"],
                        "direction": o["direction"], "weight": o["weight"]} for o in objs],
        "objectives_stats": stats,
        "ranked": ranked,
        "excluded": excluded,
        "top_id": top_id,
        "verdict": _verdict(ranked, excluded, method, top_id),
    }


def _dominates(a: List[float], b: List[float]) -> bool:
    """a 支配 b? a 在每维 ≥ b 且至少一维严格 >。"""
    at_least_one = False
    for x, y in zip(a, b):
        if x < y:
            return False
        if x > y:
            at_least_one = True
    return at_least_one


def _verdict(ranked: Optional[List[Dict[str, Any]]] = None,
             excluded: Optional[List[Dict[str, Any]]] = None,
             method: str = "", top_id: Optional[str] = None,
             reason: str = "") -> Dict[str, Any]:
    n = len(ranked) if ranked else 0
    ex = len(excluded) if excluded else 0
    narrative = (
        f"多目标排序层(方法={method}): 对 {n} 个通过硬门的候选做优先级排序"
        f"{', 剔除 ' + str(ex) + ' 个硬门失败候选' if ex else ''}。"
        f"{('排首位候选: ' + str(top_id) + '。') if top_id else '无通过候选。'}"
        f"排序仅为**实验优先级建议**(先验证哪几个), 是启发式而非活性/功能预测;"
        f"任一方向的收益结论须以湿实验 kcat/KM 与酰基酶捕获质谱收口, 不得由排序分"
        f"直接推断活性升降。{('注意: ' + reason) if reason else ''}")
    return {
        "ok": n > 0,
        "applicable": True,
        "narrative": narrative,
        "banned_phrases": [
            "排序第一 = 活性最高/功能最优",
            "加权分高 = 实验必成功",
            "帕累托前沿 = 一定可行",
            "排序分可直接推断活性升降",
        ],
        "allowed_phrases": [
            "排序仅为实验优先级建议(先验证哪几个候选)",
            "最终结论须湿实验 kcat/KM 与酰基酶捕获质谱收口",
            "硬门失败候选被剔除出排序, 不参与优先级比较",
        ],
    }


__all__ = ["ranking_layer", "make_blocker_predicate"]
