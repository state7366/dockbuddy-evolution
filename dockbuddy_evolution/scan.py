# -*- coding: utf-8 -*-
"""libdock.scan —— 突变扫描编排(B类, #312)。

================================================================================
用途: 读取一个目录下「各候选突变已算好的指标 JSON」, 聚合成候选列表, 调
libdock.ranking.ranking_layer 做多目标优先级排序, 落盘 scan_ranking.json +
字段报告。

输入约定(每条候选一个 JSON, 文件名任意, 建议 <突变名>.json):
  { "id": "F205G",                       # 缺则取文件名 stem
    "metrics": { "ddg": 1.2, "orient": 0.61, "cavity_delta_pct": 23.2, ... },  # 主路径
    "blockers": { "structure_valid": true, "consensus_hit": 0 },  # 硬门(可选) }

若 JSON 没有顶层 metrics 子 dict, scan 会按已知形态自动抽取:
  - foldx_ddg_summary(单突变): 抽 ddg=ddG_kcal, foldx_ok=ok
  - 其他: 退化为「顶层数值标量」即指标。

诚实声明: scan 只是把分散的分析结果**汇总排序**, 不重新计算任何量; 排序的
科学含义完全由喂入的 metrics 与其来源 verdict 决定。verdict 强制"优先级 ≠ 结论"。
================================================================================
"""
from __future__ import annotations

import glob
import json
import os
from typing import Any, Dict, List, Optional

import logging as _logging
LOG = _logging.getLogger("dockbuddy_evolution")
from .ranking import make_blocker_predicate, ranking_layer


def _extract_metrics(d: Dict[str, Any]) -> Dict[str, float]:
    """从单候选 JSON 抽取数值指标 dict。优先 metrics 子 dict。"""
    if isinstance(d.get("metrics"), dict):
        return {k: v for k, v in d["metrics"].items()
                if isinstance(v, (int, float)) and not isinstance(v, bool)}
    # foldx_ddg_summary 单突变: 抽 ddg + ok 标志
    if d.get("command") == "foldx-ddg" and d.get("mutations"):
        ms = d["mutations"]
        if len(ms) == 1:
            m = ms[0]
            return {"ddg": m.get("ddG_kcal"),
                    "foldx_ok": int(bool(m.get("ok")))}
    # 退化: 顶层数值标量(跳过元字段)
    out: Dict[str, float] = {}
    for k, v in d.items():
        if k in ("id", "metrics", "verdict", "note", "command", "blockers",
                 "meta", "generated_at", "ok", "available"):
            continue
        if isinstance(v, (int, float)) and not isinstance(v, bool):
            out[k] = float(v)
    return out


def scan_mutants(directory: str, *, objectives: List[Dict[str, Any]],
                weights: Optional[Dict[str, float]] = None,
                method: str = "both",
                exclude: Optional[Dict[str, List[str]]] = None,
                on_missing: str = "neutral",
                recursive: bool = False,
                out_json: Optional[str] = None
                ) -> Dict[str, Any]:
    """扫描目录下各候选 JSON → 多目标排序 → 落盘。

    参数:
      directory   候选 JSON 所在目录
      objectives  目标定义(透传给 ranking_layer)
      weights     可选权重覆盖
      method      weighted | pareto | both
      exclude     硬门 dict {"require_true":[...], "require_false":[...]}, 透传
                  make_blocker_predicate(命中者单列 excluded)
      on_missing  neutral | worst | exclude
      recursive   是否递归子目录
      out_json    输出路径(默认 <directory>/scan_ranking.json)

    返回 ranking_layer 的完整输出(含 ranked/excluded/verdict), 并附
    scanned_files 列表。
    """
    if not os.path.isdir(directory):
        return {"ok": False, "error": f"目录不存在: {directory}", "ranked": [],
                "excluded": [], "verdict": _verdict([], reason="目录不存在")}

    pattern = os.path.join(directory, "**", "*.json") if recursive \
        else os.path.join(directory, "*.json")
    files = sorted(glob.glob(pattern, recursive=recursive))
    if not files:
        return {"ok": False, "error": f"目录下无 JSON: {directory}",
                "ranked": [], "excluded": [], "verdict": _verdict([], reason="无候选")}

    require_true = list((exclude or {}).get("require_true") or [])
    require_false = list((exclude or {}).get("require_false") or [])
    pred = make_blocker_predicate(require_true=require_true,
                                  require_false=require_false) if (require_true or require_false) \
        else None

    items: List[Dict[str, Any]] = []
    scanned: List[str] = []
    for f in files:
        try:
            with open(f, encoding="utf-8") as fh:
                d = json.load(fh)
        except (OSError, ValueError) as e:
            LOG.warning("scan 跳过无法解析的 JSON: %s (%s)", f, e)
            continue
        if not isinstance(d, dict):
            continue
        iid = d.get("id") or os.path.splitext(os.path.basename(f))[0]
        metrics = _extract_metrics(d)
        item: Dict[str, Any] = {"id": iid, "metrics": metrics}
        # 透传硬门字段(顶层或 blockers 子 dict)
        for bk in set(require_true + require_false):
            if bk in d:
                item[bk] = d[bk]
            elif isinstance(d.get("blockers"), dict) and bk in d["blockers"]:
                item[bk] = d["blockers"][bk]
        items.append(item)
        scanned.append(f)

    if not items:
        return {"ok": False, "error": "未能从任何 JSON 解析出候选", "ranked": [],
                "excluded": [], "verdict": _verdict([], reason="无有效候选")}

    result = ranking_layer(items, objectives=objectives, weights=weights,
                           method=method, exclude_predicate=pred,
                           on_missing=on_missing)
    result["scanned_files"] = scanned
    result["scan_directory"] = os.path.abspath(directory)

    out_path = out_json or os.path.join(directory, "scan_ranking.json")
    try:
        os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump({k: v for k, v in result.items() if k != "item"}, fh,
                      ensure_ascii=False, indent=2)
        # 字段报告(缺 doc_type 时静默跳过, 不影响 JSON)
        try:
            import importlib
            rr = importlib.import_module("dockbuddy_dock.run_report")
            p = rr.write_run_report(out_path, result, "scan_ranking")
            if p:
                result["report_md"] = p
        except Exception as e:  # noqa: BLE001
            LOG.warning("scan 字段报告生成失败(不影响 JSON): %s", e)
        result["out_json"] = out_path
    except OSError as e:
        LOG.warning("scan 落盘失败: %s", e)
        result["out_json_error"] = str(e)

    return result


def _verdict(ranked: List[Dict[str, Any]], reason: str = "") -> Dict[str, Any]:
    n = len(ranked)
    narrative = (
        f"突变扫描编排: 汇总候选指标做多目标排序(优先级启发式, 非活性预测)。"
        f"共 {n} 个候选参与排序。排序仅为实验优先级建议, 结论须以湿实验收口; "
        f"scan 不重新计算任何量, 排序含义完全由喂入 metrics 与其来源 verdict 决定。"
        f"{('注意: ' + reason) if reason else ''}")
    return {
        "ok": n > 0, "applicable": True, "narrative": narrative,
        "banned_phrases": [
            "scan 排序第一 = 最优突变/活性最高",
            "scan 汇总排序 = 实验结论",
        ],
        "allowed_phrases": [
            "scan 仅汇总分散分析结果做优先级排程",
            "排序结论须以湿实验 kcat/KM 收口",
        ],
    }


__all__ = ["scan_mutants"]
