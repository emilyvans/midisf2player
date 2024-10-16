#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <ios>
#include <iostream>
#include <math.h>
#include <portaudio.h>
#include <vector>

std::vector<std::complex<double>> dft(std::vector<std::complex<double>> X);
std::vector<std::complex<double>> fft(std::vector<std::complex<double>> X);

float SinWave(float freq, double deltaFrame) {
	return sin(freq * 2.0 * M_PI * deltaFrame);
}

int main(int argc, char **argv) {
	// PaError init_result = Pa_Initialize();
	// if (init_result != paNoError) {
	//		std::cout << "error: " << Pa_GetErrorText(init_result) << "\n";
	//		return init_result;
	//	}
	//	for (int i = 0; i < argc; i++) {
	//		std::cout << "arg " << i + 1 << ": " << argv[i] << "\n";
	//	}

	std::vector<std::complex<double>> sin;
	float num_sample = 48000.0f;
	sin.reserve(num_sample);

	for (int i = 0; i < num_sample; i++) {
		sin.push_back(std::complex((double)SinWave(660, i / num_sample)) +
		              std::complex((double)SinWave(220, i / num_sample)));
	}

	std::vector<std::complex<double>> dft_res = fft(sin);
	std::vector<double> results;

	double max_result = 0;

	for (std::complex<double> result_complex : dft_res) {
		double result = std::abs(result_complex);
		if (result > max_result) {
			max_result = result;
		}
		results.push_back(result);
	}

	for (int i = 0; i < results.size(); i++) {
		results[i] = results[i] / max_result;
	}

	for (int i = 0; i < results.size() / 2; i++) {
		if (results[i] < 0.5)
			continue;
		std::cout << "k:" << i << std::fixed << std::setprecision(3)
							<< " normalized magnitude: " << results[i] << "\n";
	}

	// Pa_Terminate();
}

std::vector<std::complex<double>> dft(std::vector<std::complex<double>> input) {
	int N = input.size();
	int K = N;

	std::complex<double> intSum;

	std::vector<std::complex<double>> output;

	for (int k = 0; k < K; k++) {
		intSum = std::complex(0.0, 0.0);
		for (int n = 0; n < N; n++) {
			double realPart = cos(((2 * M_PI) / N) * k * n);
			double imagPart = sin(((2 * M_PI) / N) * k * n);
			std::complex<double> w(realPart, -imagPart);
			intSum += input[n] * w;
		}
		output.push_back(intSum);
	}

	return output;
}

std::vector<std::complex<double>> fft(std::vector<std::complex<double>> input) {
	int N = input.size();
	int K = N;

	std::vector<std::complex<double>> output;
	std::vector<std::complex<double>> even;
	std::vector<std::complex<double>> odd;

	if (K <= 1) {
		return input;
	}

	output.reserve(N);
	even.reserve(N / 2);
	odd.reserve(N / 2);

	for (int i = 0; i < N; i++) {
		output.push_back(NULL);
	}

	for (int i = 0; 2 * i < N; i++) {
		even.push_back(input[2 * i]);
		odd.push_back(input[2 * i + 1]);
	}

	std::vector<std::complex<double>> even_output = fft(even);
	std::vector<std::complex<double>> odd_output = fft(odd);
	for (int k = 0; k < K / 2; k++) {
		std::complex<double> t = std::polar(1.0, 2 * M_PI * k / N) * odd_output[k];
		output[k] = even_output[k] + t;
		output[N / 2 + k] = even_output[k] - t;
	}

	return output;
}
