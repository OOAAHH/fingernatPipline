#!/usr/bin/env python3
"""
批量运行fingeRNAt.py的脚本
模仿命令: python code/fingeRNAt.py -r example_inputs/1aju_model1.pdb -l example_inputs/ligands.sdf
对Das_07.pdb和extracted_poses目录下的每个sdf文件分别进行计算
"""

import os
import subprocess
import sys
from pathlib import Path
import glob

def main():
    # 定义路径
    workspace_dir = Path("/Users/joseperezmartinez/docs/rdock")
    receptor_file = workspace_dir / "Das_07.pdb"
    extracted_poses_dir = workspace_dir / "extracted_poses"
    output_base_dir = workspace_dir / "fingernaat_results"
    
    # 检查必要文件是否存在
    if not receptor_file.exists():
        print(f"错误: 受体文件不存在: {receptor_file}")
        sys.exit(1)
    
    if not extracted_poses_dir.exists():
        print(f"错误: 配体目录不存在: {extracted_poses_dir}")
        sys.exit(1)
    
    # 创建输出目录
    output_base_dir.mkdir(exist_ok=True)
    
    # 获取所有sdf文件
    sdf_files = list(extracted_poses_dir.glob("*.sdf"))
    
    if not sdf_files:
        print(f"错误: 在 {extracted_poses_dir} 中没有找到sdf文件")
        sys.exit(1)
    
    print(f"找到 {len(sdf_files)} 个sdf文件")
    
    # 处理每个sdf文件
    successful_runs = 0
    failed_runs = 0
    
    for sdf_file in sdf_files:
        print(f"\n处理文件: {sdf_file.name}")
        
        # 创建输出子目录，保留原始文件名信息
        sdf_basename = sdf_file.stem  # 不包含扩展名的文件名
        output_subdir = output_base_dir / f"Das_07_vs_{sdf_basename}"
        output_subdir.mkdir(exist_ok=True)
        
        # 构建fingeRNAt.py命令
        # 使用绝对路径避免相对路径问题
        fingernaat_script = workspace_dir / "code" / "fingeRNAt.py"
        
        cmd = [
            "python", fingernaat_script,
            "-r", str(receptor_file),
            "-l", str(sdf_file)
        ]
        
        try:
            # 在输出目录中运行命令
            result = subprocess.run(
                cmd,
                cwd=output_subdir,
                capture_output=True,
                text=True,
                timeout=300  # 5分钟超时
            )
            
            if result.returncode == 0:
                print(f"✓ 成功处理: {sdf_file.name}")
                successful_runs += 1
                
                # 保存命令输出
                with open(output_subdir / "fingernaat_output.log", "w") as f:
                    f.write(f"Command: {' '.join(cmd)}\n")
                    f.write(f"Return code: {result.returncode}\n")
                    f.write(f"STDOUT:\n{result.stdout}\n")
                    f.write(f"STDERR:\n{result.stderr}\n")
            else:
                print(f"✗ 处理失败: {sdf_file.name}")
                print(f"错误代码: {result.returncode}")
                print(f"错误信息: {result.stderr}")
                failed_runs += 1
                
                # 保存错误信息
                with open(output_subdir / "fingernaat_error.log", "w") as f:
                    f.write(f"Command: {' '.join(cmd)}\n")
                    f.write(f"Return code: {result.returncode}\n")
                    f.write(f"STDOUT:\n{result.stdout}\n")
                    f.write(f"STDERR:\n{result.stderr}\n")
                    
        except subprocess.TimeoutExpired:
            print(f"✗ 处理超时: {sdf_file.name}")
            failed_runs += 1
        except Exception as e:
            print(f"✗ 处理异常: {sdf_file.name} - {str(e)}")
            failed_runs += 1
    
    # 总结
    print(f"\n处理完成:")
    print(f"成功: {successful_runs}")
    print(f"失败: {failed_runs}")
    print(f"总计: {len(sdf_files)}")
    print(f"结果保存在: {output_base_dir}")

if __name__ == "__main__":
    main() 