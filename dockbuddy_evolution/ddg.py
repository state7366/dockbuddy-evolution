# -*- coding: utf-8 -*-
"""B3: ΔΔG 预筛 —— 便宜的"稳定性代理", 不是真 ΔΔG (2026-10-03)
=====================================================================
动机
----
定向进化一轮要筛几十~几百个突变体。每个都跑完整对接太慢, 先用一个**便宜的预筛**
把明显 destabilizing 的突变剔掉, 是标准做法。商业工具 FoldX 干这个最常用, 但本机
没有(也不打算硬依赖)。

这里给什么(**以及为什么不是朴素的总能量差**)
------------------------------------------------
最朴素的写法是 ``E(mutant) − E(WT)``, 但**这是错的**: 势能是广延量, 把 Phe 换成
Gly 会少掉 7 个重原子, 这些原子的键合/非键/GB 贡献被整个带走。实测某
Phe→Gly 突变: 总能量差 **−57.89 kcal/mol**(Phe→Ile −95.42) —— 单残基突变不可能有
这个量级, 差值的绝大部分只是"少了几个原子", 与稳定性无关。

所以本模块比较的是**同一个原子子集**: 把突变残基从两个结构里一并剔掉,

    ΔΔG_proxy = E(突变体结构 − 突变残基) − E(WT 结构 − 突变残基)

两边原子数**完全相同**, 差值才真正是"换掉这个残基之后, 蛋白其余部分的应变与
去溶剂化变化"。总能量差仍会算出来, 但只作为 ``ddG_total_kcal`` 诊断项, 明确标注
"含原子数差异, 不可作稳定性代理"。

必须反复强调的口径(写论文时别越界)
----------------------------------
1. **这不是热力学 ΔΔG**: 没有系综平均、没有熵项、**没有 unfolded 态参考**。
   真正的 ΔΔG = G_fold(mut) − G_unfold(mut) − [G_fold(wt) − G_unfold(wt)],
   而这里只有 G_fold 一侧的单结构能量。它和实验 ΔΔG 的相关性通常只有 r≈0.4~0.6。
2. 所以它的**正当用途是排序预筛**(剔除明显 destabilizing 的候选), 不是预测数值。
   报告里写 "ΔΔG_proxy" 或 "stability proxy", 绝不写 "ΔΔG = +1.2 kcal/mol"。
3. 原子数不同的两个结构, 能量**绝对值是不同尺度的**(粒子数不同 → 总能量不同数量级?
   不, 能量是广延量, 差一个侧链会让绝对值差很多)。所以:
   * 只在**同一 WT 的不同突变体之间**比较 ΔΔG_proxy(内部可比);
   * 跨蛋白比较无意义;
   * 若 |n_atom 差| 异常大(> 单个侧链量级), 说明结构不一致, 结果作废(本模块会报)。
4. 单位: OpenMM kJ/mol, 对外按 4.184 换 kcal/mol。

⚠⚠ 2026-10-04 实测标定: 本指标在本机**不可用**(噪声地板判死刑)
------------------------------------------------------------
写完代码、过了 11 项单测、真机也出了数 —— 但**出数不等于可用**。按 A1 的同款思路
做了"协议地板标定"(找一个理想值应为 0 的对照组: **同一序列的两个 AF 模型**
m4 vs m0/m1/m2, 它们叠合后 Cα RMSD 仅 0.07~0.11 Å):

====================  ====================  ====================  ====================
对照组(理想 = 0)      minimize_steps=50     minimize_steps=0      收敛到容差(-1)
====================  ====================  ====================  ====================
m4 vs m0              +92.35                −23.03                **+15.31**
m4 vs m1              +10.97                +171.01               **−23.88**
m4 vs m2              −120.29               −303.04               (未跑)
单结构耗时            ~155 s                ~35 s                 **722 ~ 1189 s**
====================  ====================  ====================  ====================

结论(**两条独立的否证, 任一条都足以判死**):
1. **噪声地板 ±20 kcal/mol**(最好的收敛模式)仍远大于单残基突变的真实量级
   (±0~5)。实测 −11.38 / −18.62(Phe→Gly / Phe→Ile 两个变体)**完全落在噪声里**。
2. **单结构 12~20 分钟**, 一轮几十个突变体要十几个小时 —— 作为"便宜的预筛",
   成本前提本身就不成立。

为什么收敛也救不回来: 蛋白势能面是**崎岖**的, 两个略有不同的结构会收敛到**不同的
局部 Basin**。比较绝对势能缺的不只是"收敛", 还缺 unfolded 态参考与系综平均 ——
这两样都不是加迭代步数能补上的。

为什么纯单点最差(±300): 未最小化的结构里有大量 0.001 Å 量级的原子冲突
(PDB 坐标只有 3 位小数 + AF 模型间微小差异), 一次 LJ 排斥就能贡献几十 kcal/mol。
所以"干脆不最小化"是**最糟**的选项, 不是省时的选项。

**处置**: 本模块保留(代码/测试/标定记录都在), 但:
  * ``ddg_proxy`` 的默认值改为 ``minimize_steps=-1``(收敛), 因为 ``>0`` 已被证明
    产出的是随机数;
  * 结果里带 ``noise_floor_kcal`` 与 ``below_noise_floor`` —— 差值落在地板内时
    明确标记"不可解释为稳定性差异";
  * **当前不得写入任何报告/论文**。要真正做 ΔΔG 预筛, 装 FoldX 或 Rosetta
    (见下)再来, 并且**先做 1~2 个有实验 ΔΔG 的突变体做方法学标定**。

为什么不用 FoldX / Rosetta
--------------------------
本机只有 openmm 8.6.1 / openbabel 3.2.1 / pdbfixer 1.12.0(见
``libdock/ligparam.py`` 的环境说明)。FoldX 需独立许可与二进制, Rosetta/PyRosetta
需编译与许可 —— 都不在。若将来装上, 建议**替换**本模块而不是并列(两套口径混用
是常见错误)。
"""

# 协议噪声地板(kcal/mol): 见上表, 收敛最小化下"同一序列的两个结构"差值的量级。
# |ΔΔG_proxy| 低于此值时, 差值与 0 无法区分。
NOISE_FLOOR_KCAL = 20.0
import logging
import os
from typing import Any, Dict, Optional

from dockbuddy_dock.solvents import OBC_ID, attach_obc_force, obc_params

LOG = logging.getLogger("libdock.ddg")

KJ_PER_KCAL = 4.184
# 单个侧链替换的原子数变化量级; 超出说明两个结构不是"同一个蛋白换一个残基"
_MAX_ATOM_DIFF = 40


def protein_energy(pdb_path: str, *, minimize_steps: int = -1,
                   solvent: str = "obc", seed: int = 20260401) -> Optional[float]:
    """蛋白单构象的 OBC 能量(kJ/mol): 补氢 → 建体系 → 局域最小化 → 单点。

    ``minimize_steps`` 语义:
      * ``> 0``  最多迭代这么多步(**不保证收敛** —— 见下);
      * ``= 0``  纯单点, 不最小化;
      * ``< 0``  **迭代到容差收敛为止**(OpenMM ``maxIterations=0`` = 不限步数)。

    ⚠ 为什么这个参数的取值会决定结果**可用还是不可用**: 局域最小化是"从起点往
    下爬", 未收敛时能量差的绝大部分是"各自爬了多深", 不是结构差异。实测同一序列
    的两个 AF 模型(m4 vs m0/m1/m2, Cα RMSD 仅 0.02~0.05 Å)在 50 步下的
    ΔΔG_proxy 是 +92.4 / +11.0 / −120.3 kcal/mol —— 噪声地板比单残基突变的真实
    效应(±0~5)大两个数量级。**未收敛的最小化等于随机数**。要用就必须收敛。

    失败返回 None(调用方降级, 别崩整条筛选链)。
    """
    import openmm as mm
    import openmm.app as app
    if not os.path.exists(pdb_path):
        LOG.warning("ddg: 结构不存在 %s", pdb_path)
        return None
    try:
        pdb = app.PDBFile(pdb_path)
        ff = app.ForceField("amber14-all.xml")
        mod = app.Modeller(pdb.topology, pdb.positions)
        mod.addHydrogens(ff, pH=7.0)
        system = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff,
                                 constraints=app.HBonds)
        if solvent == "obc" and attach_obc_force(system, LOG) is None:
            return None
        integ = mm.LangevinIntegrator(300.0 * mm.unit.kelvin,
                                      1.0 / mm.unit.picosecond,
                                      0.002 * mm.unit.picoseconds)
        integ.setRandomNumberSeed(int(seed))
        ctx = mm.Context(system, integ, mm.Platform.getPlatformByName("CPU"))
        ctx.setPositions(mod.positions)
        if minimize_steps != 0:
            mm.LocalEnergyMinimizer.minimize(
                ctx, tolerance=10.0 * mm.unit.kilojoule_per_mole / mm.unit.nanometer,
                maxIterations=max(0, int(minimize_steps)))
        return ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
            mm.unit.kilojoule_per_mole)
    except Exception as e:                                    # pragma: no cover
        LOG.warning("ddg: 能量计算失败(%s): %s", pdb_path, e)
        return None


def _n_atoms(pdb_path: str) -> int:
    n = 0
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM", "HETATM")):
                n += 1
    return n


def _residue_atom_counts(pdb_path: str) -> Dict[str, int]:
    """``{(chain, resSeq): 原子数}`` —— 用来自动定位"哪个残基被换了"。"""
    d: Dict[str, int] = {}
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if not ln.startswith(("ATOM", "HETATM")):
                continue
            key = (ln[21].strip(), ln[22:27].strip())
            d[key] = d.get(key, 0) + 1
    return d


def find_mutated_residue(wt_pdb: str, mut_pdb: str) -> Optional[str]:
    """自动定位突变残基: 两个结构里 **原子数发生变化** 的那个 ``(chain,resSeq)``。

    找不到(两侧完全一致, 或差异不止一处)返回 None —— 此时调用方必须显式给定位。
    """
    a, b = _residue_atom_counts(wt_pdb), _residue_atom_counts(mut_pdb)
    diff = [k for k in set(a) | set(b) if a.get(k) != b.get(k)]
    return diff[0] if len(diff) == 1 else None


def write_residue_excluded_pdb(pdb_path: str, out_path: str,
                               residue_key: str) -> Optional[str]:
    """写出**剔除指定残基**后的 PDB(``residue_key`` = ``(chain,resSeq)`` 字符串)。

    剔除是为了让 WT / 突变体两侧落到**同一个原子子集**上; 残基本身的原子数不同,
    不剔除的话总能量差会被"少了几个原子"主导(见模块 docstring 的实测)。
    """
    chain, seq = residue_key
    kept = []
    with open(pdb_path, encoding="utf-8", errors="ignore") as f:
        for ln in f:
            if ln.startswith(("ATOM", "HETATM")):
                if (ln[21].strip(), ln[22:27].strip()) == (chain, seq):
                    continue
            kept.append(ln.rstrip("\n"))
    if not any(l.startswith(("ATOM", "HETATM")) for l in kept):
        return None
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(kept) + "\n")
    return out_path


def ddg_proxy(wt_pdb: str, mut_pdb: str, *, minimize_steps: int = -1,
              solvent: str = "obc", seed: int = 20260401,
              residue_key=None, workdir: str = "") -> Dict[str, Any]:
    """单构象稳定性代理: **剔除突变残基后的同原子子集**能量差。

    ``residue_key`` = ``(chain, resSeq)``; None 时自动定位(两结构间原子数发生
    变化的残基)。定位不到就报 ``ok=False``, 不会拿一个错残基去算。

    返回 dict: ``{ok, ddG_proxy_kcal, ddG_proxy_kJ, E_wt_kJ, E_mut_kJ,
    ddG_total_kcal(诊断, 含原子数差异), n_atoms_subset, residue_key, ...}``。
    """
    import tempfile
    out: Dict[str, Any] = {"ok": False, "wt_pdb": wt_pdb, "mut_pdb": mut_pdb,
                           "solvent": (OBC_ID if solvent == "obc" else "vacuum")}
    for p in (wt_pdb, mut_pdb):
        if not os.path.exists(p):
            out["note"] = f"结构不存在: {p}"
            return out

    key = residue_key or find_mutated_residue(wt_pdb, mut_pdb)
    if key is None:
        out["note"] = ("无法唯一定位突变残基(两结构间原子数变化的残基不是恰好一个), "
                       "请显式传 residue_key=(chain, resSeq)")
        return out
    out["residue_key"] = list(key) if not isinstance(key, str) else key

    wd = workdir or tempfile.mkdtemp(prefix="ddg_")
    os.makedirs(wd, exist_ok=True)
    sub_wt = write_residue_excluded_pdb(wt_pdb, os.path.join(wd, "wt_sub.pdb"), key)
    sub_mut = write_residue_excluded_pdb(mut_pdb, os.path.join(wd, "mut_sub.pdb"), key)
    if not sub_wt or not sub_mut:
        out["note"] = "剔除突变残基后没有剩余原子, 无法构造可比较的子集"
        return out
    n_wt, n_mut = _n_atoms(sub_wt), _n_atoms(sub_mut)
    out["n_atoms_subset"] = {"wt": n_wt, "mut": n_mut}
    if n_wt != n_mut:
        out["note"] = (f"剔除突变残基后两侧原子数仍不等({n_wt} vs {n_mut}): "
                       f"两个结构的其余部分不一致, 差值不可用作稳定性代理")
        return out
    if _n_atoms(wt_pdb) - n_wt <= 0:
        out["note"] = "突变残基在 WT 中不存在, 请检查 residue_key"
        return out

    e_wt = protein_energy(sub_wt, minimize_steps=minimize_steps,
                          solvent=solvent, seed=seed)
    e_mut = protein_energy(sub_mut, minimize_steps=minimize_steps,
                           solvent=solvent, seed=seed)
    if e_wt is None or e_mut is None:
        out["note"] = "能量计算失败(WT 或突变体子集), 该突变体的稳定性代理不可用"
        return out
    d = e_mut - e_wt
    out.update({"E_wt_kJ": e_wt, "E_mut_kJ": e_mut,
                "ddG_proxy_kJ": d, "ddG_proxy_kcal": d / KJ_PER_KCAL,
                "n_atom_diff_full": _n_atoms(mut_pdb) - _n_atoms(wt_pdb)})
    out["ok"] = True
    out["model"] = ("单构象 OBC(GBneck2) 能量差(剔除突变残基的同原子子集): "
                    "无系综平均、无熵项、无 unfolded 态参考; "
                    "仅供同一 WT 下多突变体排序预筛")
    out["obc_params"] = obc_params() if solvent == "obc" else {}
    out["noise_floor_kcal"] = NOISE_FLOOR_KCAL
    out["below_noise_floor"] = abs(out["ddG_proxy_kcal"]) < NOISE_FLOOR_KCAL
    out["note"] = ("ΔΔG_proxy > 0 表示剔除残基后蛋白其余部分能量更高(更 destabilizing)。"
                   "只在同一 WT 的不同突变体之间比较; 不要写成实验 ΔΔG。")
    if out["below_noise_floor"]:
        out["unreliable_reason"] = (
            f"|ΔΔG_proxy| = {abs(out['ddG_proxy_kcal']):.2f} kcal/mol 低于协议噪声地板 "
            f"{NOISE_FLOOR_KCAL:.0f} kcal/mol(同一序列两个 AF 模型的对照实测), "
            f"该差值与 0 无法区分 —— 不可用于排序, 更不可写进报告")
        LOG.warning("ddg_proxy: %s", out["unreliable_reason"])
    if minimize_steps > 0:
        out["unreliable_reason"] = (
            f"minimize_steps={minimize_steps} 未收敛: 实测同一序列两结构在此模式下差值可达 "
            f"±100 kcal/mol, 差值绝大部分是'各自爬了多深'。请用 minimize_steps=-1(收敛)")
        LOG.warning("ddg_proxy: %s", out["unreliable_reason"])
    return out


def _safe_name(name: str) -> str:
    """突变体名 → 可作目录名的安全串(A123G / A123:Y:10A 都要能落盘)。"""
    return "".join(c if (c.isalnum() or c in "-_") else "_" for c in str(name))


def rank_mutants(wt_pdb: str, mutants: Dict[str, str], **kw) -> Dict[str, Any]:
    """批量: 返回按 ΔΔG_proxy 升序(更稳定的在前)的排序与逐项结果。

    每个突变体给**独立 workdir**(剔除残基后的子集 PDB 会留在里面), 便于事后核对
    "到底剔了哪个残基" —— 稳定性代理这种近似量, 可追溯性比省磁盘重要。
    """
    items = []
    base_wd = kw.pop("workdir", "") or ""
    for name, path in mutants.items():
        wd = os.path.join(base_wd, _safe_name(name)) if base_wd else ""
        r = ddg_proxy(wt_pdb, path, workdir=wd, **kw)
        r["name"] = name
        items.append(r)
    ok = [r for r in items if r.get("ok")]
    ok.sort(key=lambda r: r["ddG_proxy_kcal"])
    out = {"wt_pdb": wt_pdb, "n": len(items), "n_ok": len(ok),
           "ranking": [r["name"] for r in ok], "items": items,
           "noise_floor_kcal": NOISE_FLOOR_KCAL,
           "note": "升序 = 更稳定(ΔΔG_proxy 更小)在前"}
    spread = (max(r["ddG_proxy_kcal"] for r in ok)
              - min(r["ddG_proxy_kcal"] for r in ok)) if len(ok) > 1 else 0.0
    if len(ok) > 1 and spread < NOISE_FLOOR_KCAL:
        out["ranking_unreliable"] = True
        out["unreliable_reason"] = (
            f"全部 {len(ok)} 个突变体的 ΔΔG_proxy 极差仅 {spread:.2f} kcal/mol, "
            f"小于协议噪声地板 {NOISE_FLOOR_KCAL:.0f} kcal/mol —— 这个排序是噪声, 不可用于筛选")
        LOG.warning("rank_mutants: %s", out["unreliable_reason"])
    return out


__all__ = ["KJ_PER_KCAL", "NOISE_FLOOR_KCAL", "protein_energy", "ddg_proxy",
           "rank_mutants"]
