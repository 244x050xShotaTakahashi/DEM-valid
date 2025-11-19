#!/bin/bash
# 安息角と摩擦係数の関係を調べる数値実験実行スクリプト

set -e  # エラーで停止

# 色付き出力
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=====================================${NC}"
echo -e "${BLUE}  安息角測定実験 - 摩擦係数スタディ  ${NC}"
echo -e "${BLUE}=====================================${NC}"
echo ""

# 作業ディレクトリの確認
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# 実行ファイルの確認
if [ ! -f "./build/dem_valid" ]; then
    echo -e "${RED}エラー: 実行ファイルが見つかりません: ./build/dem_valid${NC}"
    echo "まず 'make' を実行してプログラムをビルドしてください。"
    exit 1
fi

# 結果ディレクトリの作成
RESULTS_DIR="results/friction_study"
mkdir -p "$RESULTS_DIR"

# 摩擦係数の配列
FRICTION_COEFFS=(0.1 0.2 0.3 0.4 0.5 0.6)

# 各ケースを実行
echo -e "${GREEN}6つのケースを順次実行します...${NC}"
echo ""

for mu in "${FRICTION_COEFFS[@]}"; do
    echo -e "${YELLOW}========================================${NC}"
    echo -e "${YELLOW}摩擦係数 μ = ${mu} のシミュレーション開始${NC}"
    echo -e "${YELLOW}========================================${NC}"
    
    # 入力ファイルの確認
    INPUT_FILE="inputs/friction_angle_study/input_friction_${mu}.dat"
    if [ ! -f "$INPUT_FILE" ]; then
        echo -e "${RED}エラー: 入力ファイルが見つかりません: $INPUT_FILE${NC}"
        continue
    fi
    
    # ケースごとの結果ディレクトリ
    CASE_DIR="${RESULTS_DIR}/friction_${mu}"
    mkdir -p "$CASE_DIR"
    
    # シミュレーション実行
    echo "実行コマンド: ./build/dem_valid $INPUT_FILE"
    START_TIME=$(date +%s)
    
    if ./build/dem_valid "$INPUT_FILE"; then
        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        echo -e "${GREEN}シミュレーション完了 (所要時間: ${ELAPSED}秒)${NC}"
        
        # 出力ファイルを移動
        if [ -f "data/graph11.d" ]; then
            cp data/graph11.d "${CASE_DIR}/graph11.d"
            echo "データファイルをコピー: ${CASE_DIR}/graph11.d"
        fi
        
        if [ -f "data/graph21.d" ]; then
            cp data/graph21.d "${CASE_DIR}/graph21.d"
        fi
        
        if [ -f "data/backl.d" ]; then
            cp data/backl.d "${CASE_DIR}/backl.d"
        fi
        
        # 安息角の測定
        echo "安息角を測定中..."
        if python3 src/analyze_repose_angle.py \
            --data "${CASE_DIR}/graph11.d" \
            --output "$CASE_DIR" \
            --name "friction_${mu}"; then
            echo -e "${GREEN}安息角測定完了${NC}"
        else
            echo -e "${RED}警告: 安息角測定でエラーが発生しました${NC}"
        fi
        
        # アニメーション作成（オプション）
        if [ -f "src/animate_pem.py" ]; then
            echo "アニメーションを作成中..."
            if python3 src/animate_pem.py "${CASE_DIR}/graph11.d" "${CASE_DIR}/animation_friction_${mu}.gif"; then
                echo -e "${GREEN}アニメーション作成完了: ${CASE_DIR}/animation_friction_${mu}.gif${NC}"
            else
                echo -e "${YELLOW}注意: アニメーション作成でエラーが発生しました（続行します）${NC}"
            fi
        fi
        
    else
        echo -e "${RED}エラー: シミュレーションが失敗しました${NC}"
        continue
    fi
    
    echo ""
done

# 全結果を統合したグラフを作成
echo -e "${YELLOW}========================================${NC}"
echo -e "${YELLOW}結果の統合とグラフ作成${NC}"
echo -e "${YELLOW}========================================${NC}"

if python3 src/plot_friction_vs_angle.py --input "$RESULTS_DIR" --output "${RESULTS_DIR}/friction_study_summary.png"; then
    echo -e "${GREEN}統合グラフ作成完了: ${RESULTS_DIR}/friction_study_summary.png${NC}"
else
    echo -e "${RED}警告: 統合グラフの作成に失敗しました${NC}"
fi

# 統合CSVファイルを作成
echo "統合CSVファイルを作成中..."
CSV_FILE="${RESULTS_DIR}/all_results.csv"
echo "friction_coeff,left_angle,right_angle,average_angle,std_angle,left_r2,right_r2" > "$CSV_FILE"

for mu in "${FRICTION_COEFFS[@]}"; do
    RESULT_FILE="${RESULTS_DIR}/friction_${mu}/friction_${mu}_results.csv"
    if [ -f "$RESULT_FILE" ]; then
        # CSVから値を抽出
        left_angle=$(grep "^left_angle_deg," "$RESULT_FILE" | cut -d',' -f2)
        right_angle=$(grep "^right_angle_deg," "$RESULT_FILE" | cut -d',' -f2)
        avg_angle=$(grep "^average_angle_deg," "$RESULT_FILE" | cut -d',' -f2)
        std_angle=$(grep "^std_angle_deg," "$RESULT_FILE" | cut -d',' -f2)
        left_r2=$(grep "^left_r_squared," "$RESULT_FILE" | cut -d',' -f2)
        right_r2=$(grep "^right_r_squared," "$RESULT_FILE" | cut -d',' -f2)
        
        echo "${mu},${left_angle},${right_angle},${avg_angle},${std_angle},${left_r2},${right_r2}" >> "$CSV_FILE"
    fi
done

echo -e "${GREEN}統合CSVファイル作成完了: $CSV_FILE${NC}"

# 結果のサマリーを表示
echo ""
echo -e "${BLUE}=====================================${NC}"
echo -e "${BLUE}       実験完了 - 結果サマリー        ${NC}"
echo -e "${BLUE}=====================================${NC}"
echo ""

if [ -f "$CSV_FILE" ]; then
    echo "摩擦係数 | 平均安息角 [度]"
    echo "---------|------------------"
    for mu in "${FRICTION_COEFFS[@]}"; do
        avg_angle=$(grep "^${mu}," "$CSV_FILE" | cut -d',' -f4)
        if [ ! -z "$avg_angle" ]; then
            printf "  %.1f    |     %.2f\n" "$mu" "$avg_angle"
        fi
    done
fi

echo ""
echo -e "${GREEN}全ての実験が完了しました！${NC}"
echo "結果ディレクトリ: $RESULTS_DIR"
echo "統合グラフ: ${RESULTS_DIR}/friction_study_summary.png"
echo "統合データ: ${RESULTS_DIR}/all_results.csv"
echo ""



