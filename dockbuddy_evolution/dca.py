# -*- coding: utf-8 -*-
"""libdock.dca —— 从 A3M 多序列比对估计残基对共进化耦合(纯计算, numpy)。

================================================================================
诚实声明的优先级最高(见 verdict):
- 本模块实现的是 **互信息(MI)+ 平均乘积校正(APC)**, 是共进化的**间接/弱
  信号**, 远非 plmDCA / GREMLIN / CCMpred 这类基于统计物理或伪似然的精确推断。
- 单一家族 + 有限序列数下, MI 高度受**采样噪声**与**系统发育偏倚**主导;
  近邻残基(|i-j| 小)的 MI 高往往是接触图谱的卷积伪影而非真实远程耦合。
- 输出**只作理性设计的提示性线索**, 不得当作"残基 i,j 物理相互作用"或
  "结构接触"的证据, 更不得写进论文结论。要证成结构接触须走实验(双突变
  补偿 / 共进化引导的同源建模)或 plmDCA 级别推断 + 独立验证。
================================================================================

数学:
  MI(i,j) = Σ_{a,b} p_ij(a,b) · log( p_ij(a,b) / (p_i(a)·p_j(b)) )
  APC 校正:  MI_apc(i,j) = MI(i,j) − (⟨MI_i⟩·⟨MI_j⟩) / ⟨MI⟩
  (Dunn et al. 2008, 去除系统发育/量纲偏倚, 使远程真实耦合更突出)

A3M 解析口径(与 libdock.msa 一致): 首条=query(无插入段); 其余行中**小写
字符为相对 query 的插入段(跳过)**, 大写与 '-' 均计为一个比对列 —— 列号即
query 残基号(1-based)。本模块输出下标统一用 **1-based 残基号**(与 msa 模块
对齐, 便于后续与 pocket/residue 编号直接对应)。
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

import logging as _logging
LOG = _logging.getLogger("dockbuddy_evolution")

# 20 种标准氨基酸 + gap(共 21 态); 未知字符归并到 gap。
_AA = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(_AA)}
GAP_STATE = 20
N_STATES = 21


def _read_alignment(a3m_path: str) -> List[str]:
    """解析 A3M 为 match-state 序列列表(去除小写插入段, gap 保留为 '-')。

    返回每行一个字符串, 长度=match 列数; 长度以最短行为准(外部截断)。
    """
    seqs: List[str] = []
    with open(a3m_path, encoding="utf-8", errors="ignore") as f:
        for ln in f:
            ln = ln.rstrip("\r\n")
            if not ln or ln.startswith(">"):
                continue
            col = [ch for ch in ln if not ch.islower()]  # 跳过插入段
            seqs.append("".join(col))
    return seqs


def _alignment_state_matrix(seqs: List[str]):
    """把序列列表转成 (L, N) 整数态矩阵; L=最短列数, N=序列数。

    态编号: 0..19 = 20 种 AA, 20 = gap(含未知字符)。
    """
    if not seqs:
        return None
    L = min(len(s) for s in seqs)
    N = len(seqs)
    M = np.empty((L, N), dtype=np.int64)
    for c in range(L):
        col = seqs[c] if False else [s[c] for s in seqs]  # 显式逐行取, 避免切片歧义
        for r, ch in enumerate(col):
            if ch == "-":
                M[c, r] = GAP_STATE
            else:
                M[c, r] = AA_INDEX.get(ch, GAP_STATE)
    return M


def dca_matrix(a3m_path: str, *, pseudocount: float = 0.01,
               min_seq: Optional[int] = None) -> Optional[Dict[str, Any]]:
    """读 A3M → 返回 state 矩阵 M 与边际频率 p1(归一化, 含伪计数)。

    返回 {"L","N","M"(L×N int), "p1"(L×21 归一化)} 或 None(空/太短)。
    不再预存 L×L×21×21 联合张量(大蛋白会爆内存)——联合在 _mi_matrix 内逐对现算。
    """
    seqs = _read_alignment(a3m_path)
    if len(seqs) < 2:
        return None
    M = _alignment_state_matrix(seqs)
    if M is None:
        return None
    L, N = M.shape
    if min_seq and N < min_seq:
        LOG.warning("DCA: 序列数 %d < min_seq=%d, 结果不可信", N, min_seq)
    # 一维边际(归一化, 伪计数均匀铺在 21 态上)
    p1 = np.zeros((L, N_STATES), dtype=float)
    for i in range(L):
        vals, cnts = np.unique(M[i], return_counts=True)
        p1[i, vals] = cnts
    p1 = (p1 + pseudocount / N_STATES) / (N + pseudocount)
    return {"L": L, "N": N, "M": M, "p1": p1}


def _mi_matrix(L: int, M, p1, alpha: float) -> np.ndarray:
    """逐对计算互信息(自然对数, nats)。

    联合频率用**与边际乘积成比例的伪计数**平滑(标准 DCA 做法):
        p2(a,b) = (n_ij(a,b) + α·p_i(a)·p_j(b)) / (N + α)
    这样独立列(无共进化)的 p2 ≈ p_i·p_j → MI ≈ 0, 不会因零格伪计数被抬成统一高值
    (旧方案分母不一致会把所有无关对都抬到 ~1.5 nats, 淹没真实信号)。
    """
    mi = np.zeros((L, L), dtype=float)
    N = M.shape[1]
    for i in range(L):
        ci = M[i]                      # (N,) 列 i 的态序列
        pa = p1[i][:, None]            # (21,1)
        for j in range(i + 1, L):
            cj = M[j]                  # (N,)
            n2 = np.zeros((N_STATES, N_STATES), dtype=float)
            np.add.at(n2, (ci, cj), 1.0)
            pij = (n2 + alpha * np.outer(p1[i], p1[j])) / (N + alpha)
            pb = p1[j][None, :]         # (1,21)
            with np.errstate(divide="ignore", invalid="ignore"):
                logr = np.log(pij / (pa * pb))
            finite = np.isfinite(logr)
            mi_ij = float(np.sum(pij[finite] * logr[finite]))
            mi[i, j] = mi_ij
            mi[j, i] = mi_ij
    return mi


def _apc_correct(mi: np.ndarray) -> np.ndarray:
    """平均乘积校正(APC): 去除位置自相关的量纲/系统发育偏倚。"""
    row_mean = mi.mean(axis=1)
    total_mean = mi.mean()
    if total_mean <= 0:
        return mi.copy()
    apc = mi - np.outer(row_mean, row_mean) / total_mean
    return apc


def dca_pairs(a3m_path: str, *, top_n: int = 20, pseudocount: float = 0.01,
              min_seq: Optional[int] = None, min_separation: int = 4,
              apc: bool = True) -> Dict[str, Any]:
    """主入口: A3M → 按耦合强度排序的残基对列表(1-based 残基号)。

    参数:
      top_n          返回多少对(按 APC-MI 降序, 默认 20)
      pseudocount    Laplace 平滑强度(默认 0.01)
      min_seq       序列数下限(低于则仅告警, 不拒绝; 调用方决定)
      min_separation 近邻过滤: |i-j| < 此值 的对不参与排序(默认 4, 过滤接触图卷积伪影)
      apc           True 用 APC 校正后的 MI 排序; False 用原始 MI

    返回 dict: {ok, n_sequences, n_positions, method, pairs[], verdict}。
    pairs[i] = {"i":int,"j":int,"mi":float,"apc_mi":float,"gap_i":float,"gap_j":float}
    """
    mat = dca_matrix(a3m_path, pseudocount=pseudocount, min_seq=min_seq)
    if mat is None:
        return {"ok": False, "error": "A3M 空或序列不足(<2 条)",
                "pairs": [], "verdict": _verdict([], reason="空输入")}
    L, N = mat["L"], mat["N"]
    p1, M = mat["p1"], mat["M"]
    mi = _mi_matrix(L, M, p1, alpha=pseudocount)
    apc_mi = _apc_correct(mi)

    scored = apc_mi if apc else mi
    pairs: List[Dict[str, Any]] = []
    for i in range(L):
        for j in range(i + 1, L):
            if (j - i) < min_separation:
                continue
            pairs.append({
                "i": i + 1, "j": j + 1,            # 1-based
                "mi": round(float(mi[i, j]), 5),
                "apc_mi": round(float(apc_mi[i, j]), 5),
                "gap_i": round(float(p1[i, GAP_STATE]), 4),
                "gap_j": round(float(p1[j, GAP_STATE]), 4),
            })
    pairs.sort(key=lambda p: p["apc_mi" if apc else "mi"], reverse=True)
    top = pairs[:max(0, top_n)]

    method = ("互信息(MI) + 平均乘积校正(APC) 耦合排序"
              if apc else "互信息(MI) 耦合排序")
    return {
        "ok": True,
        "n_sequences": N,
        "n_positions": L,
        "pseudocount": pseudocount,
        "min_separation": min_separation,
        "apc": apc,
        "method": method + " (非 plmDCA, 仅弱/间接共进化信号)",
        "pairs": top,
        "verdict": _verdict(top, reason=""),
    }


def _verdict(pairs: List[Dict[str, Any]], reason: str = "") -> Dict[str, Any]:
    n = len(pairs)
    top_str = "; ".join(f"{p['i']}-{p['j']}(APC={p['apc_mi']})"
                        for p in pairs[:5]) if pairs else "无"
    narrative = (
        f"共进化耦合排序(替代方案, 非 plmDCA): 返回 {n} 对候选。"
        f"Top 候选: {top_str}。MI/APC 受采样噪声与系统发育偏倚主导, 近邻伪影已按"
        f" min_separation 过滤; 该排序**仅作理性设计提示**, 不得作为结构接触/物理"
        f"相互作用证据, 不得写入论文结论。{('注意: ' + reason) if reason else ''}")
    return {
        "ok": n > 0,
        "applicable": True,
        "narrative": narrative,
        "banned_phrases": [
            "DCA 命中即结构接触/DCA=共进化接触图",
            "APC-MI 高 = 残基 i,j 存在物理相互作用",
            "共进化耦合可直接推断酶功能增益",
        ],
        "allowed_phrases": [
            "DCA(MI/APC)是共进化的弱/间接信号, 仅作理性设计线索",
            "真正确认远程耦合须 plmDCA 级别推断 + 双突变补偿实验独立验证",
        ],
    }


__all__ = ["dca_matrix", "dca_pairs", "_read_alignment", "_alignment_state_matrix"]
