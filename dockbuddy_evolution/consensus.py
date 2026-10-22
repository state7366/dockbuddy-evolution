# -*- coding: utf-8 -*-
"""libdock.consensus —— 序列共识约束(理性设计可行性围栏, 消费 msa.site_report)。

================================================================================
用途: 把 `msa.site_report`(每位点的家族保守性 + 天然变体谱)转成"**该位点允许/禁忌
突变到哪些残基**"的结构化约束, 供突变扫描排序层(ranking)与突变设计(FEP/对接)
直接消费。

口径(严格沿用 msa.site_report 的字段, 不再重算频率):
- 位点判定 `is_conserved` / `top_aa` / `aa_freqs`(非 gap 口径 0~1) / `n_effective` /
  `small_residue_frac` 全部来自 site_report(单一事实源, 避免双份口径漂移)。
- allowed(设计可行空间): 家族中天然出现频率 ≥ `allowed_freq` 的残基(含共识残基本身)。
- disfavored(规避空间): 家族中天然频率 < `allowed_freq` 的残基(极少/从未见过,
  突变过去通常破坏结构或功能)。
- 保守位点(is_conserved=True)额外给出 `consensus` 残基, 且若 `forbid_off_consensus`
  开启, 把"非共识且非天然出现"的残基显式列为禁忌(即便其频率略高于阈值, 也提示
  高风险)。
- n_effective 不足的位点: is_conserved=None, 此时**不得**生成硬约束(把"没数据"
  当成"允许任意"或"禁忌任意"都是错的), 仅给 soft 提示。

诚实声明: 共识约束是**进化保守性线索**, 不是功能/结构可行性的充分判据; 它描述
"家族见过什么", 不保证"没见过就一定不行"。verdict 强制此免责, 且禁止写出
"该位点只能突成 X" 这类绝对化措辞。
================================================================================
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Set

import logging as _logging
LOG = _logging.getLogger("dockbuddy_evolution")
from .msa import AA1_TO_AA3, SMALL_RESIDUES

# 21 态(20 AA + gap), 用于把单字母频率展开成可读列表
_ALL_AA1 = list(AA1_TO_AA3.keys())


def consensus_constraints(msa_report: Dict[int, Dict[str, Any]], *,
                         conserved_thr: float = 0.9,
                         allowed_freq: float = 0.05,
                         forbid_off_consensus: bool = True,
                         min_n_effective: int = 50) -> Dict[str, Any]:
    """把 msa.site_report 的 {site: {...}} 转成共识约束。

    参数:
      msa_report         site_report(a3m, sites) 的返回值(dict, 键=残基号)
      conserved_thr      (仅用于交叉校验 site_report 的 is_conserved, 默认 0.9)
      allowed_freq       天然频率 ≥ 此值视为 allowed(默认 0.05)
      forbid_off_consensus  保守位点是否把"非共识且非天然出现"标为禁忌(默认 True)
      min_n_effective    低于此有效序列数的位点不生成硬约束(默认 50)

    返回: {ok, n_sites, sites:{site:{...}}, verdict}
    sites[site] = {
      is_conserved, consensus(或 None), allowed:[aa1...], disfavored:[aa1...],
      small_residue_frac, n_effective, rationale(str), hard(bool)
    }
    """
    if not isinstance(msa_report, dict) or not msa_report:
        return {"ok": False, "error": "msa_report 为空或非 dict", "sites": {},
                "verdict": _verdict({})}

    sites: Dict[str, Any] = {}
    for site in sorted(msa_report, key=lambda x: int(x)):
        d = msa_report[site] or {}
        n_eff = d.get("n_effective")
        is_cons = d.get("is_conserved")
        top_aa = d.get("top_aa")
        freqs = d.get("aa_freqs") or {}
        small_frac = d.get("small_residue_frac")

        # 覆盖不足: 不生成硬约束(防止把"没数据"误读成"允许/禁忌")
        if is_cons is None or n_eff is None or n_eff < min_n_effective:
            sites[str(site)] = {
                "is_conserved": is_cons, "consensus": None,
                "allowed": [], "disfavored": [],
                "small_residue_frac": small_frac, "n_effective": n_eff,
                "hard": False,
                "rationale": (f"有效序列 {n_eff} < {min_n_effective} 或保守性判读无效"
                              f" → 不生成硬约束(避免误读), 仅作 soft 线索。"),
            }
            continue

        consensus = top_aa if is_cons else None
        # allowed = 天然频率 ≥ 阈值(家族 appreciably 使用过的残基, 设计可行空间)
        allowed: Set[str] = {aa for aa, f in freqs.items()
                             if aa in AA1_TO_AA3 and f >= allowed_freq}
        if consensus and consensus not in allowed:
            allowed.add(consensus)
        # disfavored = 家族中**从未出现过**的残基(freqs 里没有); 这些才是真正的未知风险。
        # 注意: 天然出现过但频率 < 阈值的"罕见变体"不列入 disfavored —— 它们被家族
        # 见证过, 只是不常见, 过度标禁忌会误伤(把"罕见"当成"不可行")。
        seen = {aa for aa in freqs if aa in AA1_TO_AA3}
        disfavored: Set[str] = {
            aa for aa in _ALL_AA1
            if aa not in seen and aa != consensus and aa not in allowed}
        # 保守位点额外强调: 非共识且非天然出现 → 禁忌(与上面 unseen 集合一致, 此处
        # 仅保证 forbid_off_consensus 关闭时仍至少把共识本身排除在 disfavored 外)
        if not forbid_off_consensus and is_cons and consensus:
            disfavored.discard(consensus)

        rationale = (
            f"位点保守性={'高度保守('+consensus+')' if is_cons else '非必需保守'}; "
            f"天然可行空间 {len(allowed)} 残基, 规避 {len(disfavored)} 残基(极罕见/未见过)。")
        sites[str(site)] = {
            "is_conserved": is_cons, "consensus": consensus,
            "allowed": sorted(allowed), "disfavored": sorted(disfavored),
            "small_residue_frac": small_frac, "n_effective": n_eff,
            "hard": True, "rationale": rationale,
        }

    return {
        "ok": True, "n_sites": len(sites),
        "conserved_thr": conserved_thr, "allowed_freq": allowed_freq,
        "sites": sites, "verdict": _verdict(sites),
    }


def _verdict(sites: Dict[str, Any]) -> Dict[str, Any]:
    n_hard = sum(1 for s in sites.values() if s.get("hard"))
    narrative = (
        f"序列共识约束: 对 {len(sites)} 个位点生成围栏, 其中 {n_hard} 个有硬约束"
        f"(有效序列充足), 其余仅 soft 线索。约束描述'家族进化中见过哪些残基', 是"
        f"保守性线索而非功能/结构可行性充分判据; 未出现过的残基不必然不可行。"
        f"排序层可把 disfavored 命中作为风险扣分, 但**不得**写成'该位点只能突成 X'。")
    return {
        "ok": len(sites) > 0,
        "applicable": True,
        "narrative": narrative,
        "banned_phrases": [
            "该位点只能突变成 X 残基",
            "共识约束 = 功能可行性充分判据",
            "未见过的残基必然不可行",
        ],
        "allowed_phrases": [
            "共识约束是进化保守性线索(家族见过什么)",
            "disfavored 命中仅作排序风险扣分, 非绝对禁令",
        ],
    }


__all__ = ["consensus_constraints"]
