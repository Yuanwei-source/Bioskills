#!/usr/bin/env python3
"""脚本前置门禁：按依赖清单检查本步骤所需依赖，缺失则退出 3。

设计约束（避免制造新的重复实现）：探测逻辑与"每个 stage 需要什么"**只有一份**——
它们住在 `config/dependencies.json` 与 `tools/env_check.py`，本模块只做桥接，
不自己维护依赖清单，也不自己写探测代码。

用法（放在脚本真正干活之前）：

    from _deps import require_stage
    require_stage('annot_check', __file__)

缺依赖时的输出是统一的（用途 + 安装方式 + "该步骤未执行"），并以退出码 3 结束——
与 `env_check.py --stage`、`cox1_id.py` 的 3、以及缺 `jsonschema` 的 3 同义。

**不为例外开口子**：`--help` 同样需要本步骤的依赖，因为脚本模块级就可能 import
Biopython 等库；给 help 单独放行只会制造"看起来能跑、其实缺依赖"的错觉。
"""
import importlib.util
import os
import pathlib
import sys

SKILL_DIR = pathlib.Path(__file__).resolve().parent.parent
EXIT_MISSING_DEPENDENCY = 3


def _env_check():
    """加载 tools/env_check.py（唯一实现），不执行它的 main。"""
    spec = importlib.util.spec_from_file_location('mito_env_check', SKILL_DIR / 'tools' / 'env_check.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def missing_for_stage(stage):
    """[(dep_id, purpose, install)]：该 stage 当前缺失的依赖。"""
    module = _env_check()
    manifest = module.load_manifest()
    spec = module.stage_requirements(manifest, stage)
    missing = []
    for dep_id, _tier, probe_spec, purpose, install in module.dependencies_for(
            manifest, spec.get('requires') or []):
        available, _detail, _version = module.probe({'probe': probe_spec}, None, False,
                                                    want_version=False)
        if not available:
            missing.append((dep_id, purpose, install))
    return missing


def require_stage(stage, script=None, stream=None):
    """缺依赖 → 打印统一说明并 `sys.exit(3)`；就绪 → 返回。"""
    stream = stream or sys.stderr
    if os.environ.get('MITO_ENV_GATE', '').lower() in ('off', '0', 'false', 'no'):
        # 显式测试逃逸：用 stub 依赖测试脚本自身行为时使用；会告警，不静默
        print('⚠ MITO_ENV_GATE=off：已关闭依赖门禁（仅用于测试/自检）', file=stream)
        return
    try:
        missing = missing_for_stage(stage)
    except Exception as exc:  # 清单缺失/损坏属于环境问题，同样以 3 结束
        print('✗ 无法读取依赖清单（stage: %s）: %s' % (stage, exc), file=stream)
        print('  → 该步骤未执行。完整环境盘点: python3 tools/env_check.py --setup', file=stream)
        sys.exit(EXIT_MISSING_DEPENDENCY)
    if not missing:
        return
    name = pathlib.Path(script).name if script else ''
    print('✗ 本步骤缺少依赖，未执行（stage: %s%s）' % (stage, '，%s' % name if name else ''), file=stream)
    for dep_id, purpose, install in missing:
        print('    ✗ %s —— %s' % (dep_id, purpose), file=stream)
        print('      安装: %s' % install, file=stream)
    print('  → 该步骤未执行；其余不依赖这些工具的步骤不受影响。', file=stream)
    print('  → 完整环境盘点: python3 tools/env_check.py --setup', file=stream)
    sys.exit(EXIT_MISSING_DEPENDENCY)
