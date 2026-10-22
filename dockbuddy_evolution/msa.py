# -*- coding: utf-8 -*-
"""libdock.msa —— 家族 MSA 保守性统计(通用模块, 2026-10-04)
================================================================================
来源动机: 突变体分析中临时写的 a3m 列统计被证明是突变位点选择的通用前置检查
—— 每个候选突变都要回答"这个位点家族保守吗、天然存在哪些变体"(低保守 =
低设计风险)。从 tools/oneoff/ 一次性实现收编为单一事实源。

输入格式: A3M(AlphaFold Server / ColabFold 通用)
  首条序列 = query(无插入段、无 gap); 其余行中**小写字符为相对 query 的插入**
  (跳过), 大写与 '-' 均计为一个比对列 —— 列号即 query 残基号(1-based)。

典型用法(突变位点设计前):
    from libdock.msa import msa_conservation
    rep = msa_conservation("xxx_unpaired_msa.a3m", columns=(205, 204))
    # 204 位(催化 His, 应高度保守)可作统计方法的对照锚; 205 位看变异谱。

口径注意:
  n_total = 覆盖该列的序列数(含 gap); n_effective = 非 gap 序列数;
  频率一律以 n_effective 计(UniRef 片段截断会造成高 gap 率, 见返回的 gap_frac)。
"""
from collections import Counter
from typing import Dict, Iterable, List

import logging as _logging
LOG = _logging.getLogger("dockbuddy_evolution")

# 小残基口径的**全仓唯一来源**(2026-10-05 从 tools/msa_site_report.py +
# tools/oneoff/family_stats.py 的双份定义收编; 两个 tools 脚本此后 import 这里)。
SMALL_RESIDUES = {"GLY", "ALA", "SER", "CYS", "THR", "PRO"}
AA1_TO_AA3 = {"A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
              "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
              "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
              "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"}


def msa_column_composition(a3m_path: str,
                           columns: Iterable[int]) -> Dict[int, Dict]:
    """统计 A3M 指定列(= query 残基号, 1-based)的氨基酸组成。

    返回 {col: {"n_total": int, "n_effective": int, "counter": {aa: count}}}。
    列未被任何序列覆盖时 counter 为空(调用方据此判无效)。
    """
    cols: List[int] = sorted({int(c) for c in columns})
    if not cols:
        return {}
    lo, hi = min(cols), max(cols)
    counters = {c: Counter() for c in cols}
    n_seq = 0
    with open(a3m_path, encoding="utf-8", errors="ignore") as f:
        for ln in f:
            # universal newlines 已把 \r\n→\n; rstrip("\r\n") 仅防御 lone-\r
            # 老式 Mac 换行(避免 \r 被当成一个比对列计入计数器)。
            ln = ln.rstrip("\r\n")
            if not ln or ln.startswith(">"):
                continue
            n_seq += 1                      # 含 query 本身(第 1 条)
            col = 0
            for ch in ln:
                if ch.islower():
                    continue                # 插入段: 不占比对列
                col += 1
                if col > hi:
                    break
                if col >= lo and col in counters:
                    counters[col][ch.upper()] += 1
    out: Dict[int, Dict] = {}
    for c in cols:
        cnt = counters[c]
        out[c] = {"n_total": sum(cnt.values()),
                  "n_effective": sum(v for k, v in cnt.items() if k != "-"),
                  "counter": dict(cnt)}
    LOG.info("MSA 列统计: %s (共 %d 条序列)", cols, n_seq)
    return out


def msa_conservation(a3m_path: str, columns: Iterable[int],
                     conserved_thr: float = 0.9,
                     min_effective: int = 50) -> Dict[int, Dict]:
    """在 msa_column_composition 之上做保守性判读。

    返回 {col: {
        "n_total", "n_effective", "gap_frac",
        "top_aa", "top_freq",              # 以非 gap 序列计
        "is_conserved",                    # top_freq >= conserved_thr
        "aa_freqs",                        # {aa: freq}(非 gap 口径)
        "verdict",                         # 人读结论一句话
    }}。统计有效性不足(非 gap 序列 < min_effective)时 is_conserved=None 并在
    verdict 说明 —— 避免把"没数据"读成"不保守"。
    """
    comp = msa_column_composition(a3m_path, columns)
    out: Dict[int, Dict] = {}
    for c, d in comp.items():
        n_eff = d["n_effective"]
        if n_eff < min_effective:
            out[c] = {"n_total": d["n_total"], "n_effective": n_eff,
                      "gap_frac": None, "top_aa": None, "top_freq": None,
                      "is_conserved": None, "aa_freqs": {},
                      "verdict": f"非 gap 序列仅 {n_eff}(<{min_effective}), "
                                 f"保守性判读无效"}
            continue
        cnt = {k: v for k, v in d["counter"].items() if k != "-"}
        top_aa, top_n = (max(cnt.items(), key=lambda kv: kv[1])
                         if cnt else (None, 0))
        freqs = {k: round(v / n_eff, 3) for k, v in cnt.items()}
        is_cons = (top_aa is not None and top_n / n_eff >= conserved_thr)
        if top_aa is None:
            verdict = "无覆盖"
        elif is_cons:
            verdict = f"高度保守({top_aa} {top_n / n_eff:.1%})"
        else:
            verdict = (f"非必需保守(最常见 {top_aa} 仅 {top_n / n_eff:.1%}, "
                       f"变异谱 {sorted(cnt)})")
        out[c] = {
            "n_total": d["n_total"], "n_effective": n_eff,
            "gap_frac": round(d["n_total"] and
                              (d["n_total"] - n_eff) / d["n_total"] or 0, 3),
            "top_aa": top_aa, "top_freq": round(top_n / n_eff, 3),
            "is_conserved": is_cons, "aa_freqs": freqs,
            "verdict": verdict,
        }
    return out


def site_report(a3m_path: str, sites: Iterable[int],
                conserved_thr: float = 0.9,
                min_effective: int = 50) -> Dict[int, Dict]:
    """突变位点家族保守性 + 天然变体统计(2026-10-05 从 tools/msa_site_report.py 收编)。

    在 msa_conservation 之上追加两个设计判读量:
      - ``small_residue_frac``: 天然小残基(G/A/S/C/T/P)总占比 —— 缩残基突变的
        家族先例强度。**频率直接求和, 不能再除一次 n_eff**(aa_freqs 已是
        非 gap 口径 0~1; 合成 a3m 曾把 32.8% 报成 0.5%)。
      - ``top3_variants``: 最常见三个天然变体(人读用)。

    返回 {site: {…msa_conservation 全部字段…, "small_residue_frac", "top3_variants"}}。
    sites 为残基号(1-based, = query 列号)可迭代; 无覆盖列原样透传其判读。
    """
    cols = [int(c) for c in sites]
    cons = msa_conservation(a3m_path, cols, conserved_thr=conserved_thr,
                            min_effective=min_effective)
    report: Dict[int, Dict] = {}
    for c in sorted(cons):
        d = cons[c]
        freqs = d.get("aa_freqs") or {}
        # 覆盖不足(is_conserved=None)时 small_residue_frac 必须为 None ——
        # tools 时代此分支 n_eff 非 0 会算出 0.0, 把"判读无效"报成"0% 小残基"
        # (与 32.8%→0.5% 同族的口径错误, 收编时修正并有测试锁)。
        if d.get("is_conserved") is None:
            small_frac = None
        else:
            small_frac = round(sum(v for k, v in freqs.items()
                                   if AA1_TO_AA3.get(k) in SMALL_RESIDUES), 3)
        top3 = sorted(freqs.items(), key=lambda kv: -kv[1])[:3]
        row = dict(d)
        row["small_residue_frac"] = small_frac
        row["top3_variants"] = {k: round(v, 3) for k, v in top3}
        report[c] = row
    return report


__all__ = ["AA1_TO_AA3", "SMALL_RESIDUES",
           "msa_column_composition", "msa_conservation", "site_report"]
