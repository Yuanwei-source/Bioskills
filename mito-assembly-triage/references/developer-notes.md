# 开发者笔记（人工维护）

这个文件替代了已删除的自动化"代码↔文档契约审计"（`tools/audit_code_docs.py` + 基线 + 契约测试，约
460 行）。那套东西把文档检查升格成一等公民：改两行代码要同步改文档、更新基线，成本高于它防止的
问题。现在**接口、退出码、不兼容变更**记在这里，由代码评审 + 测试兜底。

判断标准：**测试和真实数据能发现的漂移才值得自动化**；文献、措辞、语义一致性靠人读。

## 1. 项目规则：什么样的改动才值得做

排名第一的标准，也是这个 skill 存在的理由：

> **一个改动必须改变"面对一个真实线粒体时，这个 skill 说什么"，或者修掉一个真实踩到的 bug。**

| 类型 | 例子 | 处置 |
|---|---|---|
| **1. 科学正确性** | partial CDS 判定、`/transl_except`、BLAST HSP 方向/边界、环化证据、基因顺序 | 优先做；必须有失败复现 + 回归 |
| **2. 用户工作流** | 安装提示、输入检查、报告可读性、报错信息 | 可做；以"用户下一步能做什么"验收 |
| **3. 维护** | 重构、测试整理、文档 | 谨慎；必须证明**降低**复杂度（净减行数/概念） |
| **4. 元基础设施** | 经验学习系统、自动晋升、审计框架、评分体系 | **默认拒绝**，除非真实需求已经出现 |

配套约束：不新增"自造术语"进用户可见输出（severity/tier/class 之类）；报告用
`发现 → 证据 → 结论 → 下一步` 的口语结构，不用内部代号。

## 2. 测试的两档节奏

```bash
# 内循环：只跑决定生物学结论的测试（~20s）
python3 -m unittest tests.test_annot_semantics tests.test_annot_cox1_edges \
                       tests.test_annotation_policy tests.test_evidence_contracts \
                       tests.test_cox1_hsp_bounds

# 提交前 / CI：全量（~65s）
python3 -m unittest discover -s tests -p 'test_*.py'
```

`experimental/`（已停用的经验子系统）**不参与默认运行**，需要时手动：
`python3 -m unittest discover -s experimental/experience/tests -p 'test_*.py'`。

真实数据验收（QJXH、CMMC）在仓库外的 `mito-realdata-validation/`，不进公共库。

## 3. 退出码

| 场景 | 码 | 含义 |
|---|---|---|
| `check_env.sh` / `env_check.py --setup\|--daily\|--stage` | `2` | 必需依赖缺失（未执行该步骤） |
| 步骤脚本缺该步依赖（`scripts/_deps.py require_stage`） | `3` | 该步骤未执行；打印用途与安装方式 |
| `annot_check.py` | `1` | 存在 ERROR 级注释问题（结论是"有问题"，不是故障） |
| `cox1_id.py` | `0`/`1`/`2`/`3` | 有结论 / 证据不足或无命中 / 被拒绝 / 网络或**结果格式故障** |
| `reference_registry.py` | `0`/`1` | 成功 / 参数或清单错误 |
| `experimental/experience/experience.py` | `0`/`1`/`3` | 通过 / 记录非法 / 缺 `jsonschema`（校验器不可用） |

`MITO_ENV_GATE=off` 只用于测试/自检，文档中禁止用于诊断。

## 4. 不兼容变更记录

| 日期 | 变更 | 影响 |
|---|---|---|
| 2026-09 | 经验/案例子系统移出主路径（`tools/experience.py` → `experimental/experience/`），`case_records` 阶段从依赖清单移除 | 这些命令不再属于默认流程；`tests/` 不再含其测试 |
| 2026-09 | 删除代码↔文档契约审计（`tools/audit_code_docs.py` 等） | 文档一致性改人工维护（本文件） |
| 2026-09 | `case-validate` 分离 `FORMAT` / `BUSINESS` 两层输出，`VALID` 不再隐含科学结论 | 判读词单独成行，`INVALID` 附 `cause:` |
| 2026-09 | `cox1_id.py`：HSP 坐标越界/非有限/方向矛盾 = **格式故障（退出 3）**，不再降级为"无命中" | 依赖退出码的脚本需区分 `1` 与 `3` |
| 2026-09 | `jsonschema` 成为必需依赖，删除了手写替代实现 | 缺库即退出 3（提示 `pip install jsonschema`） |

## 5. 目录分工

```
scripts/            诊断与修复脚本（真正的领域逻辑）
references/         领域知识：诊断树、注释质量、证据标准、参考资料策略
schemas/            记录格式定义（case / lesson / reference-registry）
tools/              仅 env_check.py（环境就绪检查）
experimental/       已停用的元系统（经验/案例记账），默认不运行
tests/              决定"面对真实线粒体是否正确"的测试
```
