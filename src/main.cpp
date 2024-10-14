#include <cmath>
#include <complex>
#include <iostream>
#include <math.h>
#include <portaudio.h>
#include <vector>

std::vector<std::complex<double>> dft(std::vector<std::complex<double>> X);

float SinWave(float freq, double deltaFrame) {
	return sin(freq * 2.0 * M_PI * deltaFrame);
}

int main(int argc, char **argv) {
	PaError init_result = Pa_Initialize();
	if (init_result != paNoError) {
		std::cout << "error: " << Pa_GetErrorText(init_result) << "\n";
		return init_result;
	}
	for (int i = 0; i < argc; i++) {
		std::cout << "arg " << i + 1 << ": " << argv[i] << "\n";
	}

	std::vector<std::complex<double>> sin;
	sin.reserve(1024);

	for (int i = 0; i < 1024; i++) {
		sin.push_back(std::complex((double)SinWave(1023, i / 1024.0)));
	}

	std::vector<std::complex<double>> dft_res = dft(sin);

	for (int i = 0; i < dft_res.size(); i++) {
		std::complex<double> m = dft_res[i];
		if (m.imag() / 1024.0 < -0.001 || m.real() / 1024.0 > 0.001)
			std::cout << "k:" << i << " real: " << m.real() / 1024.0
								<< " imag: " << m.imag() / 1024.0 << "\n";
	}

	Pa_Terminate();
}

std::vector<std::complex<double>> dft(std::vector<std::complex<double>> X) {
	int N = X.size();
	int K = N;

	std::complex<double> intSum;

	std::vector<std::complex<double>> output;
	output.reserve(K);

	for (int k = 0; k < K; k++) {
		intSum = std::complex(0.0, 0.0);
		for (int n = 0; n < N; n++) {
			double realPart = cos(((2 * M_PI) / N) * k * n);
			double imagPart = sin(((2 * M_PI) / N) * k * n);
			std::complex<double> w(realPart, -imagPart);
			intSum += X[n] * w;
		}
		output.push_back(intSum);
	}

	return output;
}
