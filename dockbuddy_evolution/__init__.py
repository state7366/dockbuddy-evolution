# -*- coding: utf-8 -*-
"""dockbuddy_evolution —— 蛋白进化/设计分析(从 DockBuddy libdock 抽出的解耦版本)。

MSA 保守性 / 序列共识约束 / 残基共进化耦合 / ΔΔG 预筛 / 热稳定性设计(FireProt) /
多目标突变候选排序 / 扫描。
"""
from .msa import msa_conservation, site_report
from .consensus import consensus_constraints
from .dca import dca_matrix
from .ddg import protein_energy, find_mutated_residue
from .fireprot import FireProtError, submit_and_wait
from .ranking import ranking_layer, make_blocker_predicate
from .scan import scan_mutants

__version__ = "0.1.0"

__all__ = ["msa_conservation", "site_report", "consensus_constraints",
           "dca_matrix", "protein_energy", "find_mutated_residue",
           "FireProtError", "submit_and_wait", "ranking_layer",
           "make_blocker_predicate", "scan_mutants", "__version__"]
