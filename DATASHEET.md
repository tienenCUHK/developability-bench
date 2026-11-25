# DevMPNN-Bench 数据集说明（Datasheet v0.1）

> ED track 要求 Croissant / Responsible-AI 字段；本文件是数据集的叙述性说明，
> 正式提交时再转成 Croissant JSON。

## 1. 动机
衡量逆折叠模型的**可开发性可控性**：模型能否按要求移动某个理化性质（表面疏水/电荷/pI/聚集），
移动的代价有多大（scTM、多样性、recovery），以及移动时**未被请求的性质漂移多少**（泄漏）。

## 2. 组成
| 组件 | 内容 | 规模 |
|---|---|---|
| 骨架集 | 35 个单体 PDB 条目（46–415 aa），chain A | 35 |
| 设计 | 10 个臂 × 35 骨架 × 8 条 | 2 688 条 |
| 折叠结果 | ESMFold v1（3 recycles）对每条设计折叠 → scTM | 2 560 次（部分臂 33 骨架） |
| 性质 | 净电荷@pH7.4、pI、GRAVY、表面 apolar、聚集 5-mer、最长 apolar、芳香族 | 每条设计 7 项 |
| 派生量 | recovery（vs 天然序列）、多样性（平均两两 Hamming） | 每 (arm,target) |

**臂**：ProteinMPNN（+ seed 1 重复）、SolubleMPNN、MembraneMPNN(global)、LigandMPNN、
以及 5 个引导臂（surf λ10、agg λ10、charge λ+3、surf10/charge−2、surf10/charge+2）。

## 3. 采集流程
1. `prep_targets.py` 从 RCSB 下载 PDB，校验 chain A 的 CA 数在 50–450；
2. `bench30_grid.py` / `agg_grid.py` 用官方权重或 `guided_run.py` 采样（T=0.1，固定 seed）；
3. `dev_fold.py` 用 ESMFold v1 折叠并算 scTM（与输入骨架的 TM-score）；
4. `analyze_bench30.py` 计算性质 / recovery / 多样性。

## 4. 许可与出处
- 骨架：RCSB PDB（公有领域，CC0）；
- 权重：ProteinMPNN / LigandMPNN（MIT）；
- 代码：本仓库待定许可证（建议 MIT）；
- 折叠：fair-esm ESMFold（MIT），含 vendored openfold（MIT）。

## 5. 已知偏差与限制
- **只有单体、可溶、实验结构**；无膜蛋白、复合物、de novo 骨架；
- **单一 oracle**（ESMFold v1），未做 Boltz-2 / Protenix-2 复核；
- 每条 (arm,target) 只有 8 条设计，统计功效有限（重复性下限 ±0.005）；
- 性质都是**计算代理**，与实测溶解度/表达的相关性未验证（文献报告 r ≤ 0.32）；
- 表面疏水的"暴露"定义基于**输入骨架**的 Cα 埋藏度（阈值 10 Å、中位数分割），不是真实 SASA。

## 6. 伦理与合规
不涉及人类受试者、动物实验或敏感数据。全部数据来自公有数据库与开源模型输出。

## 7. 维护
作者维护；issue 与更新走仓库。骨架集固定（`data/targets30.txt`），新增臂不改动已有数据。
