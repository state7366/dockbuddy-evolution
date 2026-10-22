# dockbuddy-evolution

蛋白进化/设计分析，从 DockBuddy `libdock` 抽出的解耦版本。

| 模块 | 分析 |
|---|---|
| `msa` | 家族 MSA 保守性统计 |
| `consensus` | 序列共识约束（理性设计可行性围栏） |
| `dca` | 残基对共进化耦合（DCA） |
| `ddg` | ΔΔG 预筛（稳定性代理） |
| `fireprot` | FireProt 热稳定性设计（Web 服务客户端） |
| `ranking` | 多目标突变候选排序 |
| `scan` | 突变扫描 |

## 安装

```bash
pip install dockbuddy-evolution \
  --index-url https://pypi.cnb.cool/vcinex/bioinfo/-/packages/simple
```

## 使用

```python
from dockbuddy_evolution import msa_conservation, dca_matrix, ranking_layer, submit_and_wait
```
