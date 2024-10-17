#include <bitset>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <optional>
#include <portaudio.h>
#include <stdint.h>
#include <vector>

std::vector<std::complex<double>> fft(std::vector<std::complex<double>> input);
std::vector<std::complex<double>> fft(std::vector<float> input);

uint8_t read1(std::ifstream &file) {
	char result[1];

	file.read(result, sizeof(result));

	return result[0];
}

uint16_t read2(std::ifstream &file) {
	char result[2];
	file.read(result, sizeof(result));

	uint8_t byte1 = result[0];
	uint8_t byte2 = result[1];

	return byte2 << 8 | byte1;
}

uint32_t read3(std::ifstream &file) {
	char result[3];
	file.read(result, sizeof(result));

	uint8_t byte1 = result[0];
	uint8_t byte2 = result[1];
	uint8_t byte3 = result[2];

	return byte3 << 16 | byte2 << 8 | byte3;
}

uint32_t read4(std::ifstream &file) {
	char result[4];
	file.read(result, sizeof(result));

	uint8_t byte1 = result[0];
	uint8_t byte2 = result[1];
	uint8_t byte3 = result[2];
	uint8_t byte4 = result[3];

	return byte4 << 24 | byte3 << 16 | byte2 << 8 | byte1;
}

std::vector<uint8_t> readn_little_endian(std::ifstream &file, int n) {
	char *result = new char[n];
	file.read(result, n);

	std::vector<uint8_t> result_swap;
	result_swap.resize(n);

	for (int i = n - 1, i_swap = 0; i >= 0;) {
		result_swap[i_swap] = result[i];
		i--;
		i_swap++;
	}

	return result_swap;
}

struct WAVE_FILE {
	char magic_bytes_riff[4]; // "riff"
	char magic_bytes[4];      // "WAVE"
	uint32_t file_length;
	uint16_t format_tag;        // format category
	uint16_t channels;          // number of channels
	uint32_t samples_per_sec;   // sampling rate
	uint32_t avg_bytes_per_sec; // for buffer estimation
	uint16_t block_align;       // data block size
	uint16_t bits_per_sample;   // the bits per sample
	uint32_t data_length;
	std::vector<float> samples;
};

std::optional<WAVE_FILE> read_wave_file(const char *file_path) {
	std::ifstream file(file_path, std::ios::binary);
	if (!file.good()) {
		return std::nullopt;
	}
	WAVE_FILE wave;
	file.read(wave.magic_bytes_riff, 4);
	wave.file_length = read4(file);
	file.read(wave.magic_bytes, 4);
	read4(file); // "fmt "
	read4(file); // fmt chunk sizes
	wave.format_tag = read2(file);
	if (wave.format_tag != 1) {
		std::cout << "format isn't 1 but " << wave.format_tag << "\n";
		return std::nullopt;
	}
	wave.channels = read2(file);
	wave.samples_per_sec = read4(file);
	wave.avg_bytes_per_sec = read4(file);
	wave.block_align = read2(file);
	wave.bits_per_sample = read2(file);

	uint chunk_type = read4(file);
	while (chunk_type != 0x61746164) { // "data" but reversed
		uint32_t length = read4(file);
		file.seekg((uint32_t)file.tellg() + length);
		chunk_type = read4(file);
	}
	wave.data_length = read4(file);
	int sample_increment = ceilf(wave.bits_per_sample / (float)8);
	int64_t max = 255;
	int64_t min = 0;
	if (wave.bits_per_sample > 8) {
		int range = pow(2, wave.bits_per_sample);
		max = range / 2 - 1;
		min = -1 * (range / 2);
	}
	std::cout << "min: " << min << " max: " << max << "\n";

	uint64_t bitmask = ~(uint64_t)(pow(2, wave.bits_per_sample) - 1);

	for (int i = 0; i < wave.data_length; i += sample_increment) {
		std::vector<uint8_t> raw_sample_bytes =
				readn_little_endian(file, sample_increment);
		uint64_t raw_sample = 0;
		int64_t sample_int = 0;
		float sample = 0;
		for (int i = 0; i < sample_increment; i++) {
			raw_sample = raw_sample << 8 | raw_sample_bytes[i];
		}

		if (wave.bits_per_sample > 8) {
			uint64_t sign_bit = (uint64_t)0b1 << (sample_increment * 8 - 1);

			if ((raw_sample & sign_bit) >= 1) {
				sample_int = raw_sample | bitmask;
			} else {
				sample_int = raw_sample;
			}
		} else {
			sample_int = raw_sample;
		}

		sample = (float)sample_int / (float)max;
		// sample = sample_int;

		wave.samples.push_back(sample);

		//		std::cout << "bytes: " << std::hex << (uint16_t *)raw_sample
		//							<< ", raw(hex): " << (uint64_t)raw_sample
		//							<< ", raw(dec): " << std::dec << (uint64_t)raw_sample
		//							<< ", act: " << sample << "\n";
	}

	return std::make_optional(wave);
}

uint32_t progress = 0;

std::vector<float> samples;

int paCallback(const void *inputBuffer, void *outputBuffer,
               unsigned long framesPerBuffer,
               const PaStreamCallbackTimeInfo *timeInfo,
               PaStreamCallbackFlags statusFlags, void *userData) {
	WAVE_FILE *wave = (WAVE_FILE *)userData;
	float *out = (float *)outputBuffer;

	for (int i = 0; i < framesPerBuffer; i++) {
		*out++ = samples[progress] * 0.5;
		*out++ = samples[progress] * 0.5;
		// std::cout << wave->samples[progress] << "\n";
		progress += 1;
		if (progress >= samples.size()) {
			progress = 0;
		}
	}

	return paContinue;
}

std::vector<float> fftrange(std::vector<float> values, int start, int end) {

	std::vector<float> new_vec;

	if (values.size() < end) {
		std::cout << "error: went over max index";
		return new_vec;
	}

	for (int i = start; i <= end; i++) {
		new_vec.push_back(values[i]);
	}

	std::vector<std::complex<double>> fft_results = fft(new_vec);

	std::vector<double> results;

	double max_result = 0;

	for (std::complex<double> result_complex : fft_results) {
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
		std::cout << i << ";" << results[i] << "\n";
	}
	return new_vec;
}

std::vector<float> resampleAudio(const std::vector<float> &input,
                                 double sourceRate, double targetRate) {
	double resampleFactor = targetRate / sourceRate;
	int newLength = static_cast<int>(input.size() * resampleFactor);
	std::vector<float> output(newLength);

	// Perform linear interpolation
	for (int i = 0; i < newLength; ++i) {
		// Find the corresponding input sample position (floating point)
		double inputPos = i / resampleFactor;
		int inputIndex = static_cast<int>(inputPos);

		if (inputIndex >= input.size() - 1) {
			// Handle edge case where inputPos is near the end of the input
			output[i] = input.back();
		} else {
			// Linear interpolation between two nearest samples
			double fraction = inputPos - inputIndex;
			// a+(b-a)*t
			// (1-t) * a + t * b
			output[i] =
					(1 - fraction) * input[inputIndex] + fraction * input[inputIndex + 1];
		}
	}

	return output;
}

int main() {
	std::optional<WAVE_FILE> wave_opt = read_wave_file("test.wav");
	if (!wave_opt.has_value()) {
		std::cout
				<< "coudn't parse wave file. it either doesn't exist or isn't format 1";
		std::cin.get();
		return -1;
	}
	WAVE_FILE wave = wave_opt.value();

	std::cout << "sample rate: " << wave.samples_per_sec << "\n";

	uint64_t sample_count =
			((wave.samples.size() - 1) / wave.samples_per_sec) * 48000;

	uint64_t increment = (wave.samples.size() - 1) * 0.1;

	for (int i = 0; i < wave.samples.size() - 1; i += increment + 1) {
		std::cout << "samples " << i << " to " << i + increment << "\n";
		fftrange(wave.samples, i, i + increment);
	}

	samples = resampleAudio(wave.samples, wave.samples_per_sec, 48000);

	// return 0;

	Pa_Initialize();
	PaStream *stream;
	PaStreamParameters outputParameters;

	outputParameters.device = Pa_GetDefaultOutputDevice();
	if (outputParameters.device == paNoDevice) {
		return -1;
	}

	const PaDeviceInfo *pInfo = Pa_GetDeviceInfo(Pa_GetDefaultOutputDevice());
	if (pInfo != 0) {
		printf("Output device name: '%s' samplerate: %f\n", pInfo->name,
		       pInfo->defaultSampleRate);
	}

	outputParameters.channelCount = 2;         /* stereo output */
	outputParameters.sampleFormat = paFloat32; /* 32 bit floating point output */
	outputParameters.suggestedLatency =
			Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
	outputParameters.hostApiSpecificStreamInfo = NULL;

	PaError err = Pa_OpenStream(&stream, NULL, /* no input */
	                            &outputParameters, 48000, 0,
	                            0, /* we won't output out of range samples
	                                          so don't bother clipping them */
	                            paCallback, (void *)&wave);

	if (err != paNoError) {
		Pa_CloseStream(stream);
		stream = 0;

		std::cout << "\nearly close: " << Pa_GetErrorText(err) << "\n";
		return -1;
	}

	PaError errs = Pa_StartStream(stream);
	if (errs == paNoError) {
		std::cout << "playing\n";
		Pa_Sleep(10000);
		Pa_StopStream(stream);
	}
	Pa_CloseStream(stream);
	Pa_Terminate();
}

std::vector<std::complex<double>> fft(std::vector<float> input) {
	std::vector<std::complex<double>> vector;
	for (float value : input) {
		vector.push_back(std::complex(value));
	}
	return fft(vector);
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
