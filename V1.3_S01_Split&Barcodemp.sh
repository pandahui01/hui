#!/bin/bash
set -eo pipefail
export LC_ALL=C
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/conda/lib:/usr/local/lib


# ----------------------
# 配置参数
# ----------------------
declare -A CONFIG=(
    [input_r1]="/data/input/Files/RawData/M1ML150001777_L01_read_1.fq.gz" #read1路径，记得前面加上/data  运行前请核对路径是否有问题 ！！！
    [input_r2]="/data/input/Files/RawData/M1ML150001777_L01_read_2.fq.gz" #read2路径，记得前面加上/data  运行前请核对路径是否有问题 ！！！
    [h5_files]="/data/input/Files/RawData/zhanxiaojuan/h5/A05970A2.barcodeToPos.h5 /data/input/Files/RawData/zhanxiaojuan/h5/A05970E3.barcodeToPos.h5 /data/input/Files/RawData/zhanxiaojuan/h5/D06044G2.barcodeToPos.h5 /data/input/Files/RawData/zhanxiaojuan/h5/D06054F4.barcodeToPos.h5" # 多个H5空格分隔，记得前面加上/data  运行前请核对路径是否有问题 ！！！
    [output_dir]="/data/work/results"
    [split_parts]=20           # 分片数
    [threads]=$(nproc)         # 自动获取CPU核心数
)

# 提取测序芯片号（从input_r1文件名中获取）
flow_id=$(basename "${CONFIG[input_r1]}" | sed -E 's/(.+)_read_1.fq.gz/\1/')
echo "测序芯片号为: $flow_id"

# 创建输出目录
mkdir -p "${CONFIG[output_dir]}"

# ----------------------
# 预处理：正确分割文件
# ----------------------
# 为R1和R2创建独立分片目录
split_dir_r1="${CONFIG[output_dir]}/split_r1"
split_dir_r2="${CONFIG[output_dir]}/split_r2"
rc_split_dir="${CONFIG[output_dir]}/rc_r2"
mkdir -p "$split_dir_r1" "$split_dir_r2" "$rc_split_dir"

# 提取输入文件的基本名称
r1_basename=$(basename "${CONFIG[input_r1]}")
r2_basename=$(basename "${CONFIG[input_r2]}")

# 并行分割R1和R2（保持记录对齐）
echo "=================== 分割文件 ==================="
echo "处理时间: $(date '+%Y-%m-%d %H:%M:%S')"
seqkit split2 -p ${CONFIG[split_parts]} -j ${CONFIG[threads]} \
    -O "$split_dir_r1" "${CONFIG[input_r1]}" &
seqkit split2 -p ${CONFIG[split_parts]} -j ${CONFIG[threads]} \
    -O "$split_dir_r2" "${CONFIG[input_r2]}" &
wait
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"

# ----------------------
# 新增: 一次性生成所有R2分片的反向互补序列
# ----------------------
echo "=================== 生成反向互补序列 ==================="
echo "处理时间: $(date '+%Y-%m-%d %H:%M:%S')"
for part_num in $(seq 1 ${CONFIG[split_parts]}); do
    part_tag=$(printf "%03d" "$part_num")
    r2_split_file="${split_dir_r2}/${r2_basename%.fq.gz}.part_${part_tag}.fq.gz"
    rc_file="${rc_split_dir}/${r2_basename%.fq.gz}.part_${part_tag}_rc.fq.gz"
    echo "生成反向互补: $rc_file"
    seqtk seq -r "$r2_split_file" | pigz -c > "$rc_file" &
done
wait
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"

# ----------------------
# 核心处理函数（单分片）- 修改命名格式
# ----------------------
process_part() {
    local part_num=$1
    local h5_file=$2
    local part_tag=$(printf "%03d" "$part_num")
    local chip_id=$(basename "${h5_file}" | sed -E 's/(.+)\.barcodeToPos\.h5/\1/')
    local flow_id="$3"  # 传递测序芯片号
    local r1_base="$4"  # 传递R1文件基本名称
    local r2_base="$5"  # 传递R2文件基本名称
    local threads="$6"  # 传递线程数
    local split_parts="$7"  # 传递分片数
    local output_dir="$8"  # 传递输出目录
    
    # 创建修改后的样本特定的输出目录
    local sample_output_dir="${output_dir}/${flow_id}_${chip_id}"
    mkdir -p "$sample_output_dir"
    
    # 创建日志文件
    local log_file="${sample_output_dir}/${flow_id}_${chip_id}_part_${part_tag}.log"
    
    {
        echo "=================== 开始处理 ${flow_id}_${chip_id} 分片 $part_tag ==================="
        echo "处理时间: $(date '+%Y-%m-%d %H:%M:%S')"

        # 构建文件路径
        local r1_split_file="${split_dir_r1}/${r1_base%.fq.gz}.part_${part_tag}.fq.gz"
        # 使用预先生成的反向互补文件
        local rc_r2="${rc_split_dir}/${r2_basename%.fq.gz}.part_${part_tag}_rc.fq.gz"
        # barcodemap文件输出路径
        local barcode_output="${sample_output_dir}/${flow_id}_${chip_id}_read_2.part_${part_tag}.barcode.fq.gz"
        # 配对read1文件输出路径
        local matched_output="${sample_output_dir}/${flow_id}_${chip_id}_read_1.part_${part_tag}.matched.fq.gz"

        # 检查文件是否存在
        if [[ ! -f "$rc_r2" ]]; then
            echo "错误: 找不到反向互补文件: $rc_r2"
            return 1
        fi

        # 计算线程数，确保不会为0
        local thread_count=$(( threads / split_parts ))
        if [[ $thread_count -lt 1 ]]; then
            thread_count=108
        fi

        echo "使用线程数: $thread_count"

        # Barcode映射 - 使用新的文件名格式
        /ST_BarcodeMap/ST_BarcodeMap-0.0.1 \
            --in "$h5_file" \
            --in1 "$rc_r2" \
            --in2 "$rc_r2" \
            --out "$barcode_output" \
            --mismatch 1 \
            --thread "$thread_count" \
            --umiStart -1

        # 步骤1: 从barcode映射结果中提取ID列表
        local id_file="${sample_output_dir}/${flow_id}_${chip_id}_part_${part_tag}.ids.txt"
        pigz -dc "$barcode_output" | \
            awk 'NR % 4 == 1 {line = substr($0, 2); sub(/\|\|\|.*$/, "/1", line); print line}' > "$id_file"

        # 步骤2: 使用seqtk subseq从R1文件中提取对应序列
        echo "匹配R1分片: $r1_split_file"
        pigz -dc "$r1_split_file" | \
            seqtk subseq - "$id_file" | \
            pigz -c > "$matched_output"

        # 清理ID文件
        rm -f "$id_file"

        echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
    
    } 2>&1 | tee -a "$log_file"
}

# ----------------------
# 并行执行控制 - 修改为按H5文件依次处理
# ----------------------
export -f process_part
export split_dir_r1 split_dir_r2 rc_split_dir

echo "=================== 开始处理所有H5文件 ==================="
# 逐个处理每个h5文件
for h5_file in ${CONFIG[h5_files]}; do
    chip_id=$(basename "${h5_file}" | sed -E 's/(.+)\.barcodeToPos\.h5/\1/')
    echo "=================== 处理芯片: $chip_id ==================="
    echo "开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
    
    # 创建样本主输出目录 - 新的命名格式
    sample_output_dir="${CONFIG[output_dir]}/${flow_id}_${chip_id}"
    mkdir -p "$sample_output_dir"
    
    # 对当前h5文件的所有分片进行并行处理
    parallel --verbose -j ${CONFIG[threads]} \
        --joblog "${sample_output_dir}/${flow_id}_${chip_id}_parallel.log" \
        "process_part {1} \"$h5_file\" \"$flow_id\" \"$r1_basename\" \"$r2_basename\" ${CONFIG[threads]} ${CONFIG[split_parts]} \"${CONFIG[output_dir]}\"" ::: \
        $(seq 1 ${CONFIG[split_parts]})
    
    echo "芯片 ${flow_id}_${chip_id} 的所有分片处理完成"

    # # 合并结果到样本目录中的最终文件
    # echo "合并样本 $sample_id 的结果"
    # find "${sample_output_dir}" -name "${sample_id}_part_*.matched.fq.gz" | \
    #     xargs cat > "${sample_output_dir}/${sample_id}_final.fq.gz"
    
    echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
done

echo "=================== 所有H5文件处理完成 ==================="
echo "结束时间: $(date '+%Y-%m-%d %H:%M:%S')"