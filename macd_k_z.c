#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define N 16       // 簽名波形長度（16 個晶片）
#define MAX_ITER 10 // 迭代次數
#define alpha 0.1   // 根升餘弦濾波器的滾降係數
#define Tc 1.0     // 晶片週期
#define p 0.35     // 馬可夫機率
#define oversampling (2*N)     // 采樣率（每個晶片的分段數）

#define PI 3.141592653589793

// 半正弦波形 psi(t)
double psi(double t) {
    if (fabs(t) < Tc / 2.0) {
		return cos(PI * t / Tc);
    }
    return 0.0;
}
// 初始化簽名波形 c_k^(0)(t)
void initialize_signature_waveform(double complex *k_z, double *s_k, int user_index, int length) {
    int extended_length = 2 * length; // 2N
    double step = Tc / (oversampling); // 確保為浮點運算
    int total_samples = N * oversampling;          // 更動( extended_length * oversampling )
	//printf("Tc :%f , oversampling : %d , step : %f ",Tc,oversampling,step);
	// 初始狀態：隨機選擇 +1 或 -1
    s_k[0] = (rand() % 2 ? 1 : -1);
    // 生成馬可夫序列
    for (int i = 1; i < extended_length; i++) {
        double rand_prob = (double)rand() / RAND_MAX; // 隨機機率 [0, 1)
        if (s_k[i - 1] == 1) {
            s_k[i] = (rand_prob < p) ? -1 : 1;   // +1 -> -1 的機率 p
        } else {
            s_k[i] = (rand_prob < p) ? 1 : -1;   // -1 -> +1 的機率 p
        }
    }
	// c_k^(0)(t)
    for (int i = 0; i <= total_samples; i++) {
        double t = i * step; // 當前時間
        double complex value = 0.0+0.0*I;

        // 根據公式 (6) 計算波形c_k
        for (int n = 0; n < extended_length; n++) {
            double phase_shift = pow(-1, user_index) * n * PI / 2.0; // (-1)^k n 的相位偏移
            double complex phase = cexp(I * phase_shift); // j^((-1)^k n)
			double psi_value = t - n * Tc / 2.0 ;
            //value += phase * s_k[n] * psi_value;
			value = phase * s_k[n] * psi(psi_value);
            printf("%f \n ",psi(psi_value));
			//printf(" t: %f , nTc/2 : %f , theda : %f , psi : %f \n", t , n * Tc / 2.0 , psi_value , psi(psi_value) );//查核phase內所有相關數值
			//printf(" phase : %.1f %.1fj ,s_k[%d] : %.1f , psi: %.3f \n",creal(phase),cimag(phase),n,s_k[n],psi_value);
        }

        k_z[i] = value; // 保存波形
    }
    /*for (int i=0;i<total_samples;i++)
        printf("%f %fi\n",creal(k_z[i]),cimag(k_z[i])); */
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// 振幅正規化
void normalize(double complex *k_z, int length) {
	double max ;
	double k_z_amp[length] ;
	double k_z_arg[length];
    for (int i = 0; i < length; i++) {
        k_z_amp[i] = sqrt(creal(k_z[i]) * creal(k_z[i]) + cimag(k_z[i]) * cimag(k_z[i]));     //存取 c_k振幅
		k_z_arg[i] = atan2( creal(k_z[i]) , cimag(k_z[i] ));     //存取 c_k相位
    }
	max = k_z_amp[0];
	for (int i = 0; i < length; i++){
		if(k_z_amp[i] >= max )
			max = k_z_amp[i];
	}
	for (int i = 0; i < length; i++)
		k_z_amp[i] = k_z_amp[i] / max;
	for (int i = 0; i < length; i++)
		k_z[i] = k_z_amp[i] * cexp(I * k_z_arg[i]);
}
// 根升余弦濾波器生成
void generate_rrc_filter(double *filter, int length ) {
	double Ts = N*Tc ;
    for (int i = 0; i < length; i++) {
        double t = (i - length / 2) * Tc; // 中心化時間
        if (fabs(t) < 1e-6) {
            filter[i] = 1.0; // 避免 t = 0 的分母為零
        } else if (fabs(fabs(t) - 1.0 / (2.0 * alpha)) < 1e-6) {
            filter[i] = (alpha / sqrt(2)) * ((1 + 2 / PI) * sin(PI / (4 * alpha)) +
                                             (1 - 2 / PI) * cos(PI / (4 * alpha)));
        } else {
            filter[i] = (sin(PI * t * (1 - alpha)) +
                         4 * alpha * t * cos(PI * t * (1 + alpha))) /
                        (PI * t * (1 - pow(4 * alpha * t, 2)));
        }
    }
}
// 函數：計算累積相位
void calculate_unwrapped_phase(const double real[], const double imag[], double unwrapped_phase[], int length) {
    double prev_phase = 0.0; // 初始相位
    double cumulative_phase = 0.0;

    for (int i = 0; i < length; i++) {
        // 計算當前點的瞬時相位
        double current_phase = atan2(imag[i], real[i]);

        if (i > 0) {
            // 計算相鄰點的相位差
            double phase_diff = current_phase - prev_phase;

            // 修正相位跳變
            if (phase_diff > PI) {
                phase_diff -= 2.0 * PI;
            } else if (phase_diff < -PI) {
                phase_diff += 2.0 * PI;
            }

            // 更新累積相位
            cumulative_phase += phase_diff;
        } else {
            // 對於第一點，累積相位等於瞬時相位
            cumulative_phase = current_phase;
        }

        // 存儲累積相位
        unwrapped_phase[i] = cumulative_phase;
        prev_phase = current_phase;
    }
}

// 主流程，並保存特定迭代次數的數據
void generate_ce_waveform(double complex *k_z, int length, int save_iters[], int num_saves, const char *output_file ) {
    double filter[length];
    generate_rrc_filter(filter, length ); // 生成根升余弦濾波器

    FILE *file = fopen(output_file, "w");
    if (!file) {
        perror("Failed to open output file");
        exit(EXIT_FAILURE);
    }
    fprintf(file, "Iteration,Time,Phase\n"); // 寫入表頭

    int save_index = 0;
    for (int iter = 0; iter <= MAX_ITER; iter++) {
        // 濾波（簡化的循環卷積）
        for (int i = 0; i < length; i++) {
            double complex temp = 0.0;
            for (int j = 0; j < length; j++) {
                temp += k_z[j] * filter[(i - j + length) % length];
            }
            k_z[i] = temp;
        }

        normalize(k_z, length); // 正規化幅值

        // 保存數據（僅在指定迭代次數）
        if (save_index < num_saves && iter == save_iters[save_index]) {
			double real[length];
			double imag[length];
			double unwrapped_phase[length]; // 存儲累積相位

			// 計算累積相位
			calculate_unwrapped_phase(real, imag, unwrapped_phase, length);

            printf("This is length : %d \n",length);
            for (int i = 0; i < length; i++) {
                printf("This is ( (float)(i * Tc) )/oversampling : %f \n",( (float)(i * Tc) )/oversampling);
                fprintf(file, "%d,%f,%f\n", iter, ( (float)(i * Tc) )/oversampling, unwrapped_phase[i] ) ; // 計算相位
            }
            save_index++;
        }
    }
    fclose(file);
}

int main() {
    int total_samples = N * oversampling; // 總采樣點數       // 更動( 2 * N * oversampling )
    double s_k[2 * N]; // 初始晶片序列
    double complex k_z[total_samples]; // 簽名波形
    int save_iters[] = {0, 1, 10, 100, 1000}; // 要保存的迭代次數
    int num_saves = sizeof(save_iters) / sizeof(save_iters[0]);

    srand((unsigned int)time(NULL)); // 初始化隨機數種子

    initialize_signature_waveform(k_z,s_k, 1 ,N) ;// 初始化波形（馬可夫過程）
    //generate_ce_waveform(k_z, total_samples, save_iters, num_saves,"ce_waveform_data.csv"); // 保存數據

    //printf("Data saved to ce_waveform_data.csv\n");
    return 0;
}
