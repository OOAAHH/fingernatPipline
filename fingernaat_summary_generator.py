#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
fingeRNAt结果汇总脚本
合并PZ23、das_07_predict和das_07_dock400Pose的fingeRNAt输出结果
按指定顺序排列：
1. PZ23_fingernaatoutputs (第1行)
2. das_07_predict_fingernaat_result_outputs (第2行)
3. das_07_dock400Pose_fingernaat_result_outputs (按得分排序)
"""

import os
import re
import pandas as pd
import glob
from pathlib import Path

def extract_score_from_dirname(dirname):
    """从目录名中提取得分数值"""
    # 目录名格式: Das_07_vs_pose_XXX_score_YYY
    match = re.search(r'score_([+-]?\d+\.?\d*)', dirname)
    if match:
        return float(match.group(1))
    return float('inf')  # 如果无法解析，放到最后

def extract_pose_number_from_dirname(dirname):
    """从目录名中提取pose编号"""
    # 目录名格式: Das_07_vs_pose_XXX_score_YYY
    match = re.search(r'pose_(\d+)', dirname)
    if match:
        return int(match.group(1))
    return 0

def read_fingernaat_tsv(file_path):
    """读取fingeRNAt的TSV输出文件"""
    try:
        # 尝试不同的读取方式
        # 方法1：使用默认设置
        df = pd.read_csv(file_path, sep='\t', encoding='utf-8')
        print(f"文件 {os.path.basename(file_path)} 形状: {df.shape}")
        
        # 检查数据是否正确读取
        if df.shape[0] >= 2:
            return df
        else:
            print(f"  警告: 文件只有 {df.shape[0]} 行，尝试其他方法...")
            
        # 方法2：手动读取
        with open(file_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        
        print(f"  手动读取到 {len(lines)} 行")
        
        if len(lines) >= 2:
            # 分割标题行和数据行
            header_line = lines[0].strip()
            data_line = lines[1].strip()
            
            # 使用制表符分割
            headers = header_line.split('\t')
            data = data_line.split('\t')
            
            print(f"  标题列数: {len(headers)}, 数据列数: {len(data)}")
            
            # 确保数据列数匹配
            if len(data) != len(headers):
                # 补齐或截断数据
                if len(data) < len(headers):
                    data.extend([''] * (len(headers) - len(data)))
                else:
                    data = data[:len(headers)]
            
            # 创建DataFrame
            df = pd.DataFrame([data], columns=headers)
            print(f"  手动创建DataFrame形状: {df.shape}")
            return df
        else:
            print(f"  错误: 文件行数不足")
            return None
            
    except Exception as e:
        print(f"错误：无法读取文件 {file_path}: {e}")
        try:
            # 尝试使用latin-1编码
            df = pd.read_csv(file_path, sep='\t', encoding='latin-1')
            print(f"使用latin-1编码成功读取，文件形状: {df.shape}")
            return df
        except Exception as e2:
            print(f"使用latin-1编码也失败: {e2}")
            return None

def find_tsv_files():
    """查找所有相关的TSV文件"""
    results = {
        'pz23': None,
        'das_07_predict': None,
        'dock_poses': []
    }
    
    # 1. 查找PZ23结果
    pz23_pattern = "PZ23_fingernaatoutputs/*.tsv"
    pz23_files = glob.glob(pz23_pattern)
    if pz23_files:
        results['pz23'] = pz23_files[0]
        print(f"找到PZ23结果文件: {results['pz23']}")
    else:
        print("警告: 未找到PZ23结果文件")
    
    # 2. 查找das_07_predict结果
    predict_pattern = "das_07_predict_fingernaat_result_outputs/*.tsv"
    predict_files = glob.glob(predict_pattern)
    if predict_files:
        results['das_07_predict'] = predict_files[0]
        print(f"找到das_07_predict结果文件: {results['das_07_predict']}")
    else:
        print("警告: 未找到das_07_predict结果文件")
    
    # 3. 查找dock poses结果
    dock_pattern = "das_07_dock400Pose_fingernaat_result_outputs/Das_07_vs_pose_*/outputs/*.tsv"
    dock_files = glob.glob(dock_pattern)
    
    print(f"找到 {len(dock_files)} 个dock pose结果文件")
    
    # 为每个dock文件提取信息
    for file_path in dock_files:
        # 从路径中提取目录名
        # 路径格式: das_07_dock400Pose_fingernaat_result_outputs/Das_07_vs_pose_XXX_score_YYY/outputs/...
        path_parts = Path(file_path).parts
        dirname = path_parts[-3]  # Das_07_vs_pose_XXX_score_YYY
        
        score = extract_score_from_dirname(dirname)
        pose_num = extract_pose_number_from_dirname(dirname)
        
        results['dock_poses'].append({
            'file_path': file_path,
            'dirname': dirname,
            'score': score,
            'pose_number': pose_num
        })
    
    # 按得分排序dock poses（得分越低越好，所以升序排列）
    results['dock_poses'].sort(key=lambda x: x['score'])
    
    print(f"Dock poses按得分排序:")
    for i, pose in enumerate(results['dock_poses'][:10]):  # 只显示前10个
        print(f"  {i+1}. {pose['dirname']} (得分: {pose['score']})")
    if len(results['dock_poses']) > 10:
        print(f"  ... 还有 {len(results['dock_poses']) - 10} 个pose")
    
    return results

def normalize_row_data(df, source, score=None):
    """标准化行数据格式"""
    if df is None or df.empty:
        print(f"  警告: {source} 的数据为空")
        return None
    
    # 检查文件结构
    print(f"  处理 {source}, 文件形状: {df.shape}")
    
    # 根据DataFrame的结构获取数据行
    if len(df) >= 2:
        # 标准情况：有标题行和数据行
        data_row = df.iloc[1].copy()
    elif len(df) == 1:
        # 手动创建的DataFrame情况：只有数据行
        data_row = df.iloc[0].copy()
    else:
        print(f"  错误: {source} 文件行数不足 (只有 {len(df)} 行)")
        return None
    
    # 创建新行
    normalized_row = []
    normalized_row.append(source)  # Source
    normalized_row.append(data_row.iloc[0])  # Ligand_name
    normalized_row.append(score if score is not None else 'N/A')  # Score
    
    # 添加相互作用数据（跳过Ligand_name列）
    interaction_data = data_row.iloc[1:].tolist()
    normalized_row.extend(interaction_data)
    
    print(f"  ✓ {source} 处理完成，数据长度: {len(normalized_row)}")
    return normalized_row

def create_summary_table(file_results):
    """创建汇总表格"""
    
    # 创建结果列表
    summary_rows = []
    original_columns = None  # 用于存储原始列名
    
    # 1. 处理PZ23结果（第1行）
    print("处理PZ23结果...")
    if file_results['pz23']:
        pz23_df = read_fingernaat_tsv(file_results['pz23'])
        if pz23_df is not None and not pz23_df.empty:
            # 获取原始列名（去掉Ligand_name列）
            if original_columns is None:
                original_columns = ['Source', 'Ligand_name', 'Score'] + list(pz23_df.columns[1:])
                print(f"  使用原始列名，总共 {len(original_columns)} 列")
            
        pz23_row = normalize_row_data(pz23_df, 'PZ23')
        if pz23_row:
            summary_rows.append(pz23_row)
            print("✓ 添加PZ23结果")
    
    # 2. 处理das_07_predict结果（第2行）
    print("处理das_07_predict结果...")
    if file_results['das_07_predict']:
        predict_df = read_fingernaat_tsv(file_results['das_07_predict'])
        if predict_df is not None and not predict_df.empty:
            # 如果还没有获取原始列名，从这个文件获取
            if original_columns is None:
                original_columns = ['Source', 'Ligand_name', 'Score'] + list(predict_df.columns[1:])
                print(f"  使用原始列名，总共 {len(original_columns)} 列")
                
        predict_row = normalize_row_data(predict_df, 'Das_07_Predict')
        if predict_row:
            summary_rows.append(predict_row)
            print("✓ 添加das_07_predict结果")
    
    # 3. 处理dock poses结果（按得分排序）
    print("处理dock poses结果...")
    dock_count = 0
    for i, pose_info in enumerate(file_results['dock_poses']):
        if i % 50 == 0:  # 每50个打印进度
            print(f"  处理dock pose {i+1}/{len(file_results['dock_poses'])}")
        
        dock_df = read_fingernaat_tsv(pose_info['file_path'])
        if dock_df is not None and not dock_df.empty:
            # 如果还没有获取原始列名，从这个文件获取
            if original_columns is None:
                original_columns = ['Source', 'Ligand_name', 'Score'] + list(dock_df.columns[1:])
                print(f"  使用原始列名，总共 {len(original_columns)} 列")
        
        pose_label = f"Dock_Pose_{pose_info['pose_number']}"
        dock_row = normalize_row_data(dock_df, pose_label, pose_info['score'])
        if dock_row:
            summary_rows.append(dock_row)
            dock_count += 1
    
    print(f"✓ 添加 {dock_count} 个dock pose结果")
    
    if not summary_rows:
        print("错误: 没有有效的数据行")
        return None
    
    # 检查所有行的长度是否一致
    row_lengths = [len(row) for row in summary_rows]
    if len(set(row_lengths)) > 1:
        print(f"警告: 行长度不一致: {set(row_lengths)}")
        # 使用最小长度来统一
        min_length = min(row_lengths)
        summary_rows = [row[:min_length] for row in summary_rows]
        print(f"已统一为长度: {min_length}")
    
    # 创建DataFrame
    if original_columns is None:
        print("错误: 无法获取原始列名")
        return None
    
    # 确保列名数量与数据列数一致
    num_data_cols = len(summary_rows[0])
    if len(original_columns) != num_data_cols:
        print(f"警告: 列名数量({len(original_columns)})与数据列数({num_data_cols})不一致")
        if len(original_columns) > num_data_cols:
            original_columns = original_columns[:num_data_cols]
        else:
            # 补齐列名
            for i in range(len(original_columns), num_data_cols):
                original_columns.append(f'Extra_Feature_{i-len(original_columns)+1}')
    
    print(f"最终使用 {len(original_columns)} 列")
    print(f"前5个列名: {original_columns[:5]}")
    
    summary_df = pd.DataFrame(summary_rows, columns=original_columns)
    
    return summary_df

def main():
    """主函数"""
    print("=" * 60)
    print("fingeRNAt结果汇总脚本")
    print("=" * 60)
    
    # 查找所有TSV文件
    print("\n1. 查找TSV文件...")
    file_results = find_tsv_files()
    
    # 检查是否找到了必要的文件
    found_files = 0
    if file_results['pz23']: found_files += 1
    if file_results['das_07_predict']: found_files += 1
    found_files += len(file_results['dock_poses'])
    
    if found_files == 0:
        print("错误: 未找到任何TSV文件")
        return
    
    print(f"\n总共找到 {found_files} 个有效文件")
    
    # 创建汇总表格
    print("\n2. 创建汇总表格...")
    summary_df = create_summary_table(file_results)
    
    if summary_df is None:
        print("错误: 无法创建汇总表格")
        return
    
    # 保存结果
    output_file = "fingernaat_summary_results.tsv"
    summary_df.to_csv(output_file, sep='\t', index=False)
    print(f"\n✓ 汇总结果已保存到: {output_file}")
    
    # 显示基本信息
    print(f"\n汇总表格信息:")
    print(f"  - 行数: {len(summary_df)}")
    print(f"  - 列数: {len(summary_df.columns)}")
    print(f"  - 数据来源分布:")
    
    source_counts = summary_df['Source'].value_counts()
    for source, count in source_counts.items():
        print(f"    {source}: {count}")
    
    # 显示前几行预览
    print(f"\n前5行预览:")
    preview_cols = ['Source', 'Ligand_name', 'Score']
    print(summary_df[preview_cols].head())
    
    # 显示得分分布（仅dock poses）
    dock_scores = summary_df[summary_df['Source'].str.startswith('Dock_Pose_')]['Score']
    if len(dock_scores) > 0:
        numeric_scores = pd.to_numeric(dock_scores, errors='coerce').dropna()
        if len(numeric_scores) > 0:
            print(f"\nDock poses得分统计:")
            print(f"  - 最好得分: {numeric_scores.min():.2f}")
            print(f"  - 最差得分: {numeric_scores.max():.2f}")
            print(f"  - 平均得分: {numeric_scores.mean():.2f}")
    
    print("\n" + "=" * 60)
    print("汇总完成！")
    print("=" * 60)

if __name__ == "__main__":
    main() 