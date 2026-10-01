# 样本类群与来源初筛

## 目的与边界

在正式按类群组装/注释前，核验样本数据是否与用户声明的类群相容。用户只知道“昆虫”也能开始；属种名称可补充，不作为启动前提。`--expected-taxon Insecta` 或 `昆虫` 是待核验预期，不是已鉴定身份。不要把软件的 `Arthropoda` 模式等同于昆虫鉴定。

默认使用少量本地标记发现参考＋在线 NCBI BLAST，无需部署 nt/Kraken2/SILVA 大型分类库。发现参考仅用来寻找样本序列；参考类别、名称和遗传密码表分支不能决定样本归属。样本序列上传须有任务内明确授权；授权有效时不重复询问。下载公共发现参考不上传样本。

| 入口 | 初筛证据 | 结论范围 |
|---|---|---|
| Illumina 双端 FASTQ | 分 lane 抽样、质量处理、直接恢复的 COX1；需要时可尝试核 18S | 样本抽样中标记来源相容性；不估算整份文库组成或污染比例 |
| 组装 FASTA | 每条 record 独立的蛋白同源标记发现 | 当前组装候选归属；不验证宿主或原始测序数据 |
| FASTA＋GenBank | 首先核对记录 ID 与逐碱基序列，然后独立发现标记，并对照原 CDS | 候选归属及注释标记的独立重叠支持；不认证精确边界 |

只有 GFF3 时先按 FASTA 路线独立初筛；该初筛程序不解析 GFF3 或认证其版本。后续联合注释复核仍按现有 GFF3/GenBank 交付检查执行。若同时有组装和 reads，可分别初筛并比较两份结果，不把标记差异自动归因于测序公司。

## 1. 本地准备标记

运行前按 `config/dependencies.json` 的 `identity_prepare` 或 `identity_prepare_reads` 阶段盘点环境。MitoGeneExtractor 必须是稳定版本 >=1.9.6；旧版及 beta 拒绝使用。依赖可分散在 conda 环境，设置 `CONDA_ROOT`；也可使用绝对路径 `--mge`、`--exonerate`、`--fastp`，或 `MITOGENEEXTRACTOR`、`EXONERATE`、`FASTP` 环境变量。

```bash
# FASTQ：--r1/--r2 按相同顺序重复可保留多批次/lane
CONDA_ROOT=<conda根目录> python3 scripts/identity_prepare.py \
  --r1 <R1.fastq.gz> --r2 <R2.fastq.gz> \
  --expected-taxon Insecta --sample-pairs 1000000 \
  --outdir <样本名>_mitogenome/intermediate/identity_screen/reads

# 只有 FASTA：每条记录分别处理，禁止拼接不同来源标记
CONDA_ROOT=<conda根目录> python3 scripts/identity_prepare.py \
  --fasta <组装.fasta> --expected-taxon Insecta \
  --outdir <样本名>_mitogenome/intermediate/identity_screen/assembly

# FASTA＋GenBank：来源字段不是鉴定证据
CONDA_ROOT=<conda根目录> python3 scripts/identity_prepare.py \
  --fasta <组装.fasta> --annotation <注释.gb> \
  --expected-taxon Insecta \
  --outdir <样本名>_mitogenome/intermediate/identity_screen/assembly_annotation
```

长时间命令用 `task_manager.py start` 包装；任务目录放 `intermediate/tasks/`，核对最终退出码和日志。输出目录须全新，不覆盖原数据或旧初筛。

程序对全部 lane 检查配对与 FASTQ 完整性后，固定随机种子抽样；原始文件只读。采样量影响恢复能力，100万对是工程起点，不是生物学判据。若支持不足，可在新目录扩大至500万对；达到任务预算后保持证据不足，不无限重试。全量输入核验仍需要读完整文件，不能承诺固定几分钟完成。

默认从 `config/identity_reference_sources.json` 下载固定提交和 SHA256 的10条 COX1 蛋白发现参考。下载仅发生在任务目录，不写入安装目录；可用 `--protein-ref` 提供小型广泛类群蛋白参考，但须解释其覆盖限制。参考不是所有生物的全面覆盖；默认表5与表2分支仅用于标记发现，不自动决定后续注释密码表。

FASTQ 经 fastp 处理后交给 MitoGeneExtractor。程序不用其默认共识直接判定：按原模板合并重叠 mates，mate 冲突不计支持，默认需5个模板且优势碱基比例>=0.9；显著混合列使该共识拒绝上传，并传入最终未解决事项。无支持列和 alignment 缺口不能被删掉后拼接。默认只提交连续>=300 bp 的有效片段。阈值均是可解释的工程筛选，不认证真实线粒体、排除 NUMT 或确定物种。

发现参考按 `--treat-references-as-individual` 独立运行匹配，避免同一标记的 reads 因竞争蛋白参考被拆散而造成假低支持。不同参考得到的查询及反向互补表示去重并保留来源；多个发现分支不是多个独立标记。

输出 `markers.json` 包含模式、预期类群、输入/参考/query SHA256、抽样信息、工具命令、标记来源及混合观察；`query_*.fasta` 是真正来自样本的查询片段。没有可靠标记返回1，不能推断非昆虫。

## 2. 在线 BLAST 与可恢复判读

确认任务授权后：

```bash
python3 scripts/identity_search.py \
  --markers <初筛目录>/markers.json \
  --outdir <初筛目录>/online --allow-public-upload
```

可通过 `--email` 或 `NCBI_EMAIL` 提供 NCBI 联系邮箱。只上传最多20条、每条最多5000 bp 的标记；超过预算须审阅标记，不静默截断。首轮使用不限制昆虫的 nucleotide BLAST，让非昆虫候选参与比较。

明确使用 NCBI `core_nt`，并保存 XML 中实际数据库名。NCBI 当前默认会把 `nt` 请求切换到 `core_nt`；不能把请求参数写成 nt 就宣称搜索了完整 nt。公共数据库覆盖有限，无命中不否定样本类群。

结果通过 `XML2_S` 保存原文并核对每条 query ID/长度，复用现有 COX1 HSP 解析器。相同序列对应多个 XML2 `HitDescr` 时保留全部 taxid。使用覆盖>=0.8、identity>=75%的合格命中，以及总 bitscore 距最优5%以内的命中共同计算 taxonomy 最低共同祖先（LCA）；总分仅汇总不重叠且共线的合格 HSP。此阈值用于粗筛，不能作为物种鉴定阈值，近缘参考缺失可能降低可判读等级。

若强竞争命中因重叠等歧义不能可靠汇总，不允许较弱命中单独决定方向；结果为不足。若近最佳命中达到本次30条的返回上限，也保持不足，因为未返回的同分候选可能改变 LCA。

每个搜索任务的 NCBI 请求间隔至少10秒，每个 RID 轮询至少间隔60秒。多个搜索目录没有共享请求时钟，因此在线任务须串行调度以遵守整体频率限制；本地准备可并行。RID、时钟、输入绑定和原始结果保存在 `search_state.json` 与查询目录；同目录拒绝并发查询。任务等待超出本次轮询预算时：

```bash
python3 scripts/identity_search.py \
  --markers <初筛目录>/markers.json \
  --outdir <初筛目录>/online --allow-public-upload --resume
```

resume 会核对输入、query 和已保存结果；已完成时复用绑定结果，不重复上传。若提交时断线且未取得 RID，不能自动确定服务器是否已接收；需查看 `submission.txt`/状态后处理，不盲目重新提交。服务超时、无效 XML、RID 失效或 taxonomy 不完整都不能变成“非昆虫”。

## 3. 有冲突时的可选核 18S 路线

若线粒体初筛冲突/疑似混合，需要核标记帮助判断宿主。可在新 FASTQ 准备目录增加 `--recover-18s`：只下载 Rfam RF01960 单个模型，用 Infernal 招募 reads，保留 mates，经 SPAdes 局部组装后再次用模型确认，提交相应片段。可用 `--ssu-model` 提供模型；依赖阶段为 `identity_prepare_ssu`。

这条短 reads 组合路线当前标记为实验性：模型命中只支持 RNA 家族身份，不确定宿主；未完成独立混样恢复率评估。不能以18S候选产出认证标本，未恢复也不是否定昆虫。若局部组装有歧义，保留 contig 候选及未解决状态。COX1 与另一个线粒体基因是同一来源层面的证据，不能代替核标记。

## 4. 用户结论与后续操作

`online/identity_result.json` 的总体状态：

| 状态 | 用户含义 | 后续 |
|---|---|---|
| `compatible` | 标记证据与预期类群相容 | 继续适用路线；不改变用户物种名 |
| `suspected_mismatch` | 合格标记支持不同分支，存在来源矛盾 | 扩大抽样、核标记及样本/批次核对；暂缓按原类群定稿 |
| `suspected_mixed` | 不同有效查询同时出现相容和不相容来源 | 比较各 contig/lane/标记，不能自动选择最符合预期的那个 |
| `insufficient` | 无标记、歧义命中、taxonomy/恢复证据不足 | 可带限制继续适用检查；不得判为非昆虫 |
| `service_or_format_failure` | 服务/结果故障，尚无科学判读 | 恢复或诊断服务，不能作生物学结论 |

即使 LCA 到 species，也只是公共近似参考的暂定归属，不独立确认物种。输出显式保留 `species_identified=false` 与 `company_error_established=false`。来源矛盾可能涉及样本、寄生/共生生物、混样、数据库记录或流程选择，不能自行断言公司返回错误。

报告首页写一句“样本来源初筛结论”，在“样本与数据”说明原始数据/候选范围、核与线粒体标记、参考登录号、未确认属种、冲突和未评估事项。既有 `cox1_id.py` 保持单区域在线查询能力；它要求坐标，不是 FASTQ 标记提取器。

## 依据

- [MitoGeneExtractor 论文](https://doi.org/10.1111/2041-210X.14075)；[上游实现](https://github.com/cmayer/MitoGeneExtractor)：蛋白参考对 reads/FASTA 的标记恢复，默认共识与质量/配对限制。
- [Rfam RF01960](https://rfam.org/family/RF01960)；[Infernal 使用说明](https://docs.rfam.org/en/latest/genome-annotation.html)：真核 SSU rRNA 模型；具体短 reads 组合路线仍需评估。
- [NCBI API](https://blast.ncbi.nlm.nih.gov/doc/blast-help/urlapi.html)；[使用限制](https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html)：程序化批量查询与请求间隔。
