# run_fingernaat_pipeline 使用说明

本脚本为端到端流程：从参赛提交的 PDB 中拆出模型与配体、生成配体 SDF、调用 fingeRNAt 计算相互作用指纹，并在 summarize 阶段构建 canonical 口袋索引，最终汇总为 CSV/Excel/TSV 及仅包含 PDF 的热图。

核心脚本：`bin/fingeRNAt/run_fingernaat_pipeline.py`

---

## 变更摘要（重要）

- Pipeline logic
  - 移除 CCD/XLSX 与自动分类（不再依赖 `scripts/pdb_ligand_report`）
  - 新增手工配体选择：通过 `--ligand-csv`（两列 `puzzle,label`）指定每个 puzzle 的配体标签
  - 识别规则：与任意标签匹配的残基组视为配体；未匹配 ATOM → `rna.pdb`；未匹配 HETATM 忽略（不输出 `others.pdb`）
  - 模型拆分：内置只保留前 5 个 MODEL/MODULE（不可配置）
  - 写出配体时优先使用 HETATM 行生成 `ligand.pdb`/`ligand.sdf`（若无 HETATM，退化为该组全部坐标）
  - 当 `--root` 指向 `Puzzles/original` 时，预测模型会从 `Puzzles/PZxx/step2/*.pdb` 读取，并按文件名中的 submitter 与编号组织为 `PZxx/<submitter>/model_XX`
  - `finalize` 阶段不再自动使用“第一个 solution”兜底，而是依赖 `Puzzles/solution/PZxx/split.txt` 与 `Puzzles/table/PZxx.csv` 的显式映射：每个 `(submitter, model)` 必须映射到某个 `solution_i`（0 基索引），否则记录错误并跳过对齐

- CLI options
  - 新增：`--ligand-csv`（当执行 `prep` 阶段时必填）
  - 保留：`--fingernaat-script`、`--python-exec`、`--steps`、`--overwrite`、`--dry-run`、`--no-heatmaps`
  - 移除：`--ccd`、`--ligand-guide`、`--no-others`

- Data model
  - 残基分组新增标志：`has_atom`、`has_hetatm` 用于区分聚合物/小分子来源
  - `category` 不再由 CCD 推断，写入元数据时固定为 `ligand`（仅对被标签选中的残基）
  - 取消 `suspect_missing_separator` 与所有 CCD 相关字段（如 `ccd_name`/`ccd_formula` 等）
  - `rna.pdb` 仅包含 ATOM；未匹配的 HETATM 不再输出

- How to use now
  1) 准备 CSV：两列 `puzzle,label`，每行一个标签；支持 `RES`、`RES_CHAIN`、`RES_CHAIN_RESSEQ`、`RES_CHAIN_RESSEQ_ICODE`（下划线或冒号分隔皆可）
  2) 执行预处理：`--ligand-csv` 必填，生成 `rna.pdb`、每配体 `ligand.pdb/.sdf`
  3) 执行 fingeRNAt：对每个配体调用外部 `fingeRNAt.py`
  4) 汇总结果：输出全局 CSV、每 puzzle 的交互明细 Excel、canonical 口袋 TSV、binding‑pocket TSV+Excel 以及仅 PDF 的热图

## 功能分解
- 预处理（prep）
  - 逐个 PDB 按 MODEL/MODULE 切分模型
  - 基于手工提供的配体标签 CSV（`--ligand-csv`）匹配配体；非配体的 ATOM 归入 `rna.pdb`，未匹配的 HETATM 忽略
  - 为每个模型：写出 `rna.pdb`，并为每个配体残基写出 `ligand.pdb` 与 `ligand.sdf`（需 OpenBabel/pybel）；全模型命名在最终阶段（finalize）统一处理
  - 记录每个配体的 `metadata.json`
- 运行（run）
  - 对每个 `rna.pdb` 与对应 `ligand.sdf` 调用外部 `fingeRNAt.py`（以每个配体为单位）
  - 将 fingeRNAt 的输出（tsv/csv）保存在每个配体目录的 `fingernaat/` 下，并写出 `fingernaat_status.json`
- 汇总（summarize）
  - 聚合所有配体的运行状态与输出，写出总表 `fingernaat_summary.csv`
  - 按 puzzle 生成交互明细工作簿（`…_fingernaat_results.xlsx`），含 `interactions`（交互行）与 `runs`（运行记录）两个工作表
  - 从交互明细导出 canonical 口袋坐标的标准化 TSV（`PZxx_fingernaat_summary_results.tsv`，供热图与 binding‑pocket 聚合）
  - 基于 canonical TSV 聚合 binding‑pocket 指纹，写出 `PZxx_fingernaat_summary_results_binding_pocket.tsv` 与带字符级着色的 Excel（`…_binding_pocket.xlsx`）
  - 调用内置脚本 `regenerate_interaction_summary.py` 生成仅包含 PDF 的热图（综合 / per‑interaction / per‑ligand），以及透视统计 Excel

---

## 依赖环境
- Python 3.8+
- OpenBabel 的 Python 绑定（`openbabel.pybel`），用于将配体 PDB 转 SDF；若缺失则无法产生 `ligand.sdf`，从而无法运行 fingeRNAt
- fingeRNAt 可执行脚本（例如本仓库中的 `bin/fingeRNAt/fingeRNAt.py` 或外部安装）
- 仅在 summarize 的可视化部分需要：
  - `pandas`（导出 TSV）
  - `openpyxl`（写 Excel 工作簿、binding‑pocket 指纹与透视统计）
  - `matplotlib`、`seaborn`（生成综合/分组/按相互作用拆分的 PDF 热图）

---

## 输入组织
- `--root` 指向包含若干 puzzle 子目录的根目录
  - 传统邮件提取数据：例如 `99_emailExtrctPZfiles/`，每个 puzzle 目录下直接存放原始提交的 `*.pdb`
  - 规范化预测数据：当 `--root` 为 `Puzzles/original` 时，脚本会自动从 `Puzzles/PZxx/step2/*.pdb` 读取预测模型，而不是 `Puzzles/original/PZxx`（后者仅作为原始归档）
- `--puzzles` 指定要处理的子目录名列表（如 `PZ43 PZ47 PZ49`）
- 每个 puzzle 目录下放若干 `*.pdb` 提交文件
- `--ligand-csv` 配体标签 CSV：两列 `puzzle,label`，按 puzzle 指定要作为配体的残基标签（见下文）
- step2 模型命名约定（仅当 `--root Puzzles/original` 时生效）：
  - 预测模型文件名形如 `PZxx_<submitter>_<index>.pdb`（例如 `PZ4_Adamiak_3.pdb`）
  - 预处理后写入 `output_root/PZxx/<submitter>/model_XX`，其中 `XX` 等于 `<index>`，用于和 `Puzzles/table/PZxx.csv` 中的 `(submitter, model)` 映射保持一致

---

## 运行示例

### 场景一：老版 99_emailExtrctPZfiles（仅预测）

最常见的一键运行（准备 + 运行 + 汇总）：

```bash
python bin/fingeRNAt/run_fingernaat_pipeline.py \
  --root 99_emailExtrctPZfiles \
  --puzzles PZ43 PZ47 PZ49 \
  --ligand-csv path/to/ligands.csv \
  --fingernaat-script bin/fingeRNAt/fingeRNAt.py \
  --python-exec python \
  --output-root fingernaat_pipeline_outputs
```

只做预处理（生成 `rna.pdb` 与每配体 `ligand.sdf`）：

```bash
python bin/fingeRNAt/run_fingernaat_pipeline.py \
  --root 99_emailExtrctPZfiles \
  --puzzles PZ43 \
  --ligand-csv path/to/ligands.csv \
  --fingernaat-script bin/fingeRNAt/fingeRNAt.py \
  --steps prep
```

只运行 fingeRNAt（假定预处理已完成）：

```bash
python bin/fingeRNAt/run_fingernaat_pipeline.py \
  --root 99_emailExtrctPZfiles \
  --puzzles PZ43 \
  --ligand-csv path/to/ligands.csv \
  --fingernaat-script bin/fingeRNAt/fingeRNAt.py \
  --steps run
```

只做汇总与可视化（假定 fingeRNAt 输出已存在）：

```bash
python bin/fingeRNAt/run_fingernaat_pipeline.py \
  --root 99_emailExtrctPZfiles \
  --puzzles PZ43 \
  --ligand-csv path/to/ligands.csv \
  --fingernaat-script bin/fingeRNAt/fingeRNAt.py \
  --steps summarize
```

### 场景二：规范化 step2 + solution + super 对齐（推荐）

以官方 `Puzzles/original` + `Puzzles/solution` + `Puzzles/table` 为基础，一次性完成预测/solution 预处理、fingeRNAt 运行、结果汇总与 super 对齐：

```bash
python bin/fingeRNAt/run_fingernaat_pipeline.py \
  --root Puzzles/original \
  --solutions-root Puzzles/solution \
  --puzzles PZ38 PZ43 PZ47 PZ49 PZ4 PZ21 PZ25 \
  --ligand-csv bin/fingeRNAt/firnatINPUT.csv \
  --fingernaat-script bin/fingeRNAt/fingeRNAt.py \
  --python-exec python \
  --output-root fingernaat_pipeline_outputs \
  --steps prep prep_solution run summarize finalize \
  --overwrite
```

说明：
- 预测模型将从 `Puzzles/PZxx/step2/*.pdb` 读取，并写入 `fingernaat_pipeline_outputs/PZxx/<submitter>/model_XX`
- solution 将从 `Puzzles/solution/PZxx` + `split.txt` 读取，并写入 `fingernaat_pipeline_outputs/PZxx/<solution_stem>/model_XX`
- `finalize` 阶段会按 `Puzzles/table/PZxx.csv` 与 `split.txt` 中的顺序（0 基 `solution_i`）为每个 prediction 选择对应 solution 进行 RNA super 对齐

---

## 全部参数
- `--root` 根目录（必填）。包含各 puzzle 子目录
- `--puzzles` 要处理的 puzzle 名（必填），支持多个
- `--ligand-csv` 配体标签 CSV（当执行 `prep` 阶段时必填）
- `--fingernaat-script` fingeRNAt 的可执行脚本路径（必填）
- `--python-exec` 调用 fingeRNAt 所用的 Python，可填 `python`/`python3` 或绝对路径，默认 `python`
- `--output-root` 输出根目录，默认 `fingernaat_pipeline_outputs`
- `--steps` 要执行的阶段，默认 `prep run summarize finalize`（推荐完整流程：`prep prep_solution run summarize finalize`）
- `--overwrite` 已有结果是否覆盖（默认不覆盖）
- `--dry-run` 仅打印流程，不落盘也不执行外部命令
-（已移除）CCD/自动分类与 XLSX 指南；本流程仅根据 `--ligand-csv` 匹配配体
- `--no-heatmaps` 在 summarize 阶段不生成 TSV + 热图

---

## 配体标签 CSV（--ligand-csv）使用指南

`--ligand-csv` 用于手动指定每个 puzzle 要当作配体的残基标签。支持两种格式：

1) 两列 CSV（重复行）：

```
puzzle,label
PZ43,LIG
PZ47,UNL_A_101
PZ49,LIG:B:19
```

2) 多列 CSV（同一行多个别名）：

```
PZ43,LIG,UNL,
PZ47,LIG,UNL,GNG
PZ49,LIG,UNL,2ZY
```

说明：
- 多行可为同一 puzzle 指定多个标签
- label 支持以下形式（下划线或冒号分隔均可；残基名不区分大小写）：
  - `RES`（仅按残基名匹配，如 `LIG`）
  - `RES_CHAIN` 或 `RES:CHAIN`（如 `LIG_B`）
  - `RES_CHAIN_RESSEQ` 或 `RES:CHAIN:RESSEQ`（如 `LIG_B_19`）
  - `RES_CHAIN_RESSEQ_ICODE` 或 `RES:CHAIN:RESSEQ:ICODE`
- 仅当某个残基组（在单一模型内，按 `resname,chain,resseq,icode` 聚合）与任一标签匹配时，才会被当作配体抽取
- 其他 ATOM 记录写入 `rna.pdb`；未匹配的 HETATM 记录写入 `others.pdb`

注意：
- 该模式不再使用 CCD 或自动分类；请确保标签能唯一定位目标配体
- 若标签仅给出 `RES`，且文件中存在多个同名小分子，将全部被当作配体抽取

---

## 输出目录结构（新契约）
所有产出位于 `--output-root` 下（默认 `fingernaat_pipeline_outputs/`）。在最终阶段（finalize）会重命名/对齐全模型文件：

```
output_root/
  fingernaat_summary.csv                # 全部配体级别的汇总（含运行状态、输出文件列表等）
  PZ43/
    <短名或原始PDB文件名>/
      model_01/
        rna.pdb                         # 仅聚合物（RNA/蛋白），用作对齐参考/移动体
        <stem>.pdb 或 <stem>_01.pdb     # 全模型（含配体）：
                                        #  - solution：<solution_stem>.pdb（不做对齐，仅改名）
                                        #  - prediction：<source_stem>_<model:02d>.pdb（按 PZxx.csv / split.txt 显式选择对应的 solution_i 进行 super 对齐；若该模型无有效映射，则仅改名不对齐）
        # others.pdb 不再输出（未匹配的 HETATM 被忽略）
        ligand_<RES>_<CHAIN>_<SEQ>_occN/
          ligand.pdb
          ligand.sdf                    # 需 pybel 才会生成
          metadata.json                 # 该配体的上下文元数据
          fingernaat/                   # 在此工作目录调用 fingeRNAt
            *.tsv / *.csv              # fingeRNAt 输出（不同版本位置略有不同）
            outputs/*.tsv              # 若 fingeRNAt 写入 outputs 子目录，脚本会读取这里的 TSV 做汇总
          fingernaat_status.json        # 外部命令返回码/标准输出/错误输出
    PZ43_fingernaat_results.xlsx                     # 交互行（interactions）+ 运行记录（runs）
    PZ43_fingernaat_summary_results.tsv              # 从 interactions 导出的 canonical 口袋 TSV（用于热图与 binding‑pocket）
    PZ43_fingernaat_summary_results_binding_pocket.tsv   # binding‑pocket 指纹（每残基一列，bitstring 长度 = 12）
    PZ43_fingernaat_summary_results_binding_pocket.xlsx # binding‑pocket 指纹 Excel（bit 1–12 着色，顶部内嵌 legend）
    heatmaps/                                        # 可视化输出（可选生成）
      interaction_summary_heatmap.pdf                # solution/submitter 分组的相互作用汇总热图
      comprehensive_heatmap.pdf                      # 综合热图（上：指纹；下：每残基总指纹数量曲线）
      by_interaction/*.pdf                           # 按相互作用拆分的热图
      per_ligand/*.pdf                               # 每配体的 residue×interaction 热图
      visual_pivots.xlsx                             # 汇总透视统计（per_ligand / per_residue / totals）
```

说明：
- 若 `fingernaat_status.json.returncode != 0`，请查看同目录 `stderr` 内容定位问题（缺少 SDF、fingeRNAt 依赖等）
- `summarize` 阶段仅会解析 `fingernaat/outputs/*.tsv` 内的交互行（若直接在 `fingernaat/` 下有 TSV，也会被记录路径，但不参与行级汇总）

---

## 使用要点与常见问题
- 生成 SDF 依赖 `openbabel.pybel`；若缺失将只写出 `ligand.pdb`，并在运行阶段报告 `sdf_missing`
- 重新跑同一目录请加 `--overwrite`，否则检测到已有 `fingernaat_status.json` 会跳过
- 可用 `--steps prep` 分步检查预处理输出（便于排查分类/切分问题）

---

## 参考
- 脚本入口：`bin/fingeRNAt/run_fingernaat_pipeline.py`
- 热图相关：`bin/fingeRNAt/export_heatmap_tsv.py`、`bin/fingeRNAt/regenerate_interaction_summary.py`、`bin/fingeRNAt/visualize_fingernaat_heatmap.py`

---

## Solution 管线（prep_solution）

- 目标：对 `--solutions-root/PZ{i}` 下的官方解（solution）执行与提交相同的配体抽取流程，便于后续统一汇总与对比。
- 输入组织：每个 `PZ{i}` 目录包含若干 solution PDB 以及一个 `split.txt`；`split.txt` 的每行以制表符/逗号/空格分隔，第一列为 solution PDB 的基名（如 `8ITS`，对应 `8ITS.pdb`）。
  - 注意：若按照 `split.txt` 的基名找不到对应 `<stem>.pdb`，流程将立即报错并停止（不会做模糊匹配）。
- 执行示例：
  ```bash
  python bin/fingeRNAt/run_fingernaat_pipeline.py \
    --solutions-root Puzzles/solution \
    --puzzles PZ40 PZ43 \
    --ligand-csv bin/fingeRNAt/firnatINPUT.csv \
    --fingernaat-script bin/fingeRNAt/fingeRNAt.py \
    --python-exec python \
    --output-root fingernaat_pipeline_outputs \
    --steps prep_solution
  ```
- 输出与排序：
  - 最终产出：solution 全模型为 `<solution_stem>.pdb`，prediction 全模型为 `<source_stem>_<model:02d>.pdb`；并保留 `rna.pdb` 与每配体 `ligand.pdb/.sdf`，不输出 `others.pdb`
  - 在汇总（`summarize`）阶段，solution 的记录会被优先置顶：
    - Excel 的 `interactions` 与 `runs` 工作表中，solution 行排在最前
    - 导出的按相互作用的 TSV 与 by_interaction 热图中，solution 也显示在第一行

---

## 最终阶段（finalize）说明
- 目的：
  - Solution 的 `model_full.pdb` 改名为 `<solution_stem>.pdb`（不做对齐）
  - Prediction 的 `model_full.pdb` 使用 `rna.pdb` 对齐到“该模型对应的 solution_i 的 RNA”，再将坐标变换应用到 full，输出为 `<source_stem>_<model:02d>.pdb`
  - 移除原始 `model_full.pdb`，仅保留新契约命名的全模型文件
- 参考选择与映射规则：
  - `solution_i` 由 `Puzzles/solution/PZxx/split.txt` 的行顺序定义，**i 从 0 开始**：第 1 行→`solution_0`、第 2 行→`solution_1`，依此类推
  - `Puzzles/table/PZxx.csv` 至少包含列 `0`（submitter）、`model`、`solution`，其中 `solution` 为 0 基 `solution_i`；负数表示该行不参与 super 对齐
  - 对 prediction 目录 `PZxx/<submitter>/model_XX`，finalize 使用 `(submitter, XX)` 在 `PZxx.csv` 中查找 `solution_i`，再根据 `solution_i` 选取相应 solution stem 下同一 `model_id` 的 `rna.pdb` 作为对齐参考
  - 若缺少 `split.txt`、缺少 `PZxx.csv` 或某个 `(submitter, model)` 无合法 `solution_i` / 参考 RNA，则记录错误日志并跳过该模型的对齐；**不会再自动回退到任意 solution**
- 对齐基元：主链原子（P, OP1/2/3, O5', C5', C4', C3', O3', O4', C1', C2'；兼容 * 记法），仅使用两侧共有的原子集合；若匹配点少于 3，则仅改名不对齐
- 约束：需要安装 `biopython`（`Bio.SVDSuperimposer`）与 `numpy`
