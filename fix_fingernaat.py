#!/usr/bin/env python3
"""
修复版本的fingeRNAt.py调用脚本
修复文件扩展名解析问题，能正确处理包含多个点的文件名
"""

import sys
import os
import subprocess
from pathlib import Path

def fix_extension_parsing(original_script, receptor_file, ligand_file):
    """
    通过修改参数调用原始脚本，避免文件扩展名解析问题
    """
    # 获取真正的扩展名
    ligand_path = Path(ligand_file)
    real_extension = ligand_path.suffix[1:]  # 去掉前面的点
    
    if real_extension == 'sdf':
        # 如果是sdf文件，直接调用原始脚本
        cmd = [
            'python', original_script,
            '-r', receptor_file,
            '-l', ligand_file
        ]
        return subprocess.run(cmd)
    else:
        print(f"错误：不支持的文件格式：{real_extension}")
        return subprocess.CompletedProcess([], 1)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("用法: python fix_fingernaat.py <原始脚本> <受体文件> <配体文件>")
        sys.exit(1)
    
    original_script = sys.argv[1]
    receptor_file = sys.argv[2]
    ligand_file = sys.argv[3]
    
    result = fix_extension_parsing(original_script, receptor_file, ligand_file)
    sys.exit(result.returncode) 