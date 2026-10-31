# 可克隆源码交付

本包的 `CloudDRTS_Source_v7.bundle` 是完整、可克隆的本地 Git 源码仓库，包含本次源代码、依赖锁、必要案例和历史版 ANDES 许可/源码。大体积的生成结果在外层程序包中；Git 仓库通过入口重新生成。

```bash
git clone CloudDRTS_Source_v7.bundle CloudDRTS_Source_v7
cd CloudDRTS_Source_v7
git log -1
```

先按 README 建立自编环境并运行 `CloudDRTS/run_all.py`，再建立独立 ANDES 环境运行 `external_run/run_reproduction.py`。生成阶段说明时需要两边结果及 `requirements_notes.txt`。

所有代码以相对路径定位。源码提交号记录在 `GIT_SOURCE.json`，Git bundle 完整性已经通过 `git bundle verify` 和克隆后逐文件 SHA-256 核对。

没有配置或推送在线远端。得到双方认可的 Git 远端地址后，可以在此克隆中添加 remote 并推送；本次未代为创建或公开远端仓库。
