#include "midi_note_freq.hpp"
#include <bitset>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdlib>
#include <format>
#include <fstream>
#include <iostream>
#include <optional>
#include <portaudio.h>
#include <stdint.h>
#include <string.h>
#include <vector>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

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

std::vector<uint8_t> readn_samples(std::ifstream &file, int n,
                                   int sample_size) {
	char *result = new char[n];
	file.read(result, n);

	std::vector<uint8_t> result_swap;
	result_swap.resize(n);

	for (int sample = 0; sample < n; sample += sample_size) {
		for (int i = sample_size - 1, i_swap = 0; i > 0;) {
			result_swap[sample + i_swap] = result[sample + i];
			i--;
			i_swap++;
		}
	}

	delete[] result;

	return result_swap;
}

struct WAVE_FILE_SAMPLER_CHUNK_LOOP {
	uint32_t identifier;
	uint32_t type;
	uint32_t start;
	double start_percent;
	uint32_t end;
	double end_percent;
	uint32_t fraction;
	uint32_t play_count;
};

struct WAVE_FILE_SAMPLER_CHUNK {
	uint32_t manufacturer;
	uint32_t product;
	uint32_t sample_period;
	uint32_t midi_unity_note;
	uint32_t midi_pitch_fraction;
	uint32_t smpte_format;
	uint32_t smpte_offset;
	uint32_t sample_loops_count;
	uint32_t sampler_data_length;
	std::vector<WAVE_FILE_SAMPLER_CHUNK_LOOP> Loops;
	// std::vector<uint8_t> sample_data;
};

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
	std::optional<WAVE_FILE_SAMPLER_CHUNK> sampler_chunk;
	uint32_t data_length;
	std::vector<float> samples;
};

WAVE_FILE_SAMPLER_CHUNK read_wave_file_sampler_chunk(std::ifstream &file) {
	WAVE_FILE_SAMPLER_CHUNK chunk;
	chunk.manufacturer = read4(file);
	chunk.product = read4(file);
	chunk.sample_period = read4(file);
	chunk.midi_unity_note = read4(file);
	chunk.midi_pitch_fraction = read4(file);
	chunk.smpte_format = read4(file);
	chunk.smpte_offset = read4(file);
	chunk.sample_loops_count = read4(file);
	chunk.sampler_data_length = read4(file);

	for (int i = 0; i < chunk.sample_loops_count; i++) {
		WAVE_FILE_SAMPLER_CHUNK_LOOP chunk_loop;
		chunk_loop.identifier = read4(file);
		chunk_loop.type = read4(file);
		chunk_loop.start = read4(file);
		chunk_loop.end = read4(file);
		chunk_loop.fraction = read4(file);
		chunk_loop.play_count = read4(file);
		chunk.Loops.push_back(chunk_loop);
	}

	return chunk;
}

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

	uint32_t chunk_type = read4(file);
	while (chunk_type != 0x61746164) { // "data" but reversed
		uint32_t length = read4(file);
		uint32_t start = file.tellg();
		switch (chunk_type) {
		case 0x6C706D73:
			wave.sampler_chunk = read_wave_file_sampler_chunk(file);
			break;
		}
		file.seekg(start + length);
		chunk_type = read4(file);
	}
	wave.data_length = read4(file);
	std::cout << wave.data_length << "\n";
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

	std::vector<uint8_t> samples =
	    readn_samples(file, wave.data_length, sample_increment);

	for (int i = 0; i < wave.data_length; i += sample_increment) {
		uint64_t raw_sample = 0;
		int64_t sample_int = 0;
		float sample = 0;
		for (int n = 0; n < sample_increment; n++) {
			raw_sample = raw_sample << 8 | samples[i + n];
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

		//		std::cout << "bytes: " << std::hex << (uint16_t
		//*)raw_sample
		//							<< ", raw(hex):
		//" << (uint64_t)raw_sample
		//							<< ", raw(dec):
		//" << std::dec << (uint64_t)raw_sample
		//							<< ", act: " <<
		// sample
		//<< "\n";
	}

	if (wave.sampler_chunk.has_value()) {
		WAVE_FILE_SAMPLER_CHUNK *sampler_chunk = &wave.sampler_chunk.value();
		int data_length = wave.data_length;
		for (int i = 0; i < sampler_chunk->sample_loops_count; i++) {
			sampler_chunk->Loops[i].start_percent =
			    1.0 / data_length * (float)sampler_chunk->Loops[i].start;
			sampler_chunk->Loops[i].end_percent =
			    1.0 / data_length * (float)sampler_chunk->Loops[i].end;
		}
	}

	return std::make_optional(wave);
}

uint32_t progress = 0;
double sample_rate;
bool looped = false;
double sample_offset;

std::vector<float> samples;
std::vector<float> samples_loop;

int paCallback(const void *inputBuffer, void *outputBuffer,
               unsigned long framesPerBuffer,
               const PaStreamCallbackTimeInfo *timeInfo,
               PaStreamCallbackFlags statusFlags, void *userData) {
	WAVE_FILE *wave = (WAVE_FILE *)userData;
	float *out = (float *)outputBuffer;
	if (wave->sampler_chunk == std::nullopt) {
		for (int i = 0; i < framesPerBuffer; i++) {
			*out++ = samples[progress] * 0.5;
			*out++ = samples[progress + 1] * 0.5;
			// std::cout << wave->samples[progress] << "\n";
			progress += wave->channels;
			if (progress >= samples.size()) {
				progress = 0;
			}
		}
	} else {
		WAVE_FILE_SAMPLER_CHUNK_LOOP loop =
		    wave->sampler_chunk.value().Loops[0];

		std::cout << loop.start << "(" << loop.start * loop.start_percent
		          << ")/" << loop.end << "(" << loop.end * loop.end_percent
		          << ") : " << samples.size() << ":" << wave->data_length
		          << "\n";
		for (int i = 0; i < framesPerBuffer; i++) {
			*out++ = samples[progress] * 0.5;
			*out++ = samples[progress + 1] * 0.5;
			// std::cout << wave->samples[progress] << "\n";
			progress += wave->channels;
			if ((!looped && progress >= samples.size()) ||
			    (looped && progress >= (loop.end * sample_offset))) {
				looped = true;
				progress = loop.start * sample_offset;
			}
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

	if (start > -1) {
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
                                 double sourceRate, double targetRate,
                                 int channels) {
	double resampleFactor = targetRate / sourceRate;
	int newLength =
	    static_cast<int>((input.size() / (float)channels) * resampleFactor) *
	    channels;
	std::vector<float> output(newLength);

	// Process each channel separately
	for (int ch = 0; ch < channels; ++ch) {
		for (int i = 0; i < newLength / channels; ++i) {
			// Find the corresponding input sample position
			// (floating point) for this channel
			double inputPos = i / resampleFactor;
			int inputIndex = static_cast<int>(inputPos);

			// Index of the current sample for the current channel
			int inputOffset = inputIndex * channels + ch;
			int outputOffset = i * channels + ch;

			if (inputOffset >= input.size() - channels) {
				// Handle edge case where inputPos is near the
				// end of the input for this channel
				output[outputOffset] = input[input.size() - channels + ch];
			} else {
				// Linear interpolation between two nearest
				// samples for this channel
				double fraction = inputPos - inputIndex;
				output[outputOffset] = (1 - fraction) * input[inputOffset] +
				                       fraction * input[inputOffset + channels];
			}
		}
	}

	return output;
}

int main(int argc, char **argv) {
	int sample_rate_default = 48000;
	std::optional<WAVE_FILE> wave_opt = read_wave_file("test.wav");
	if (!wave_opt.has_value()) {
		std::cout << "coudn't parse wave file. it either doesn't exist or "
		             "isn't format 1";
		return -1;
	}
	WAVE_FILE wave = wave_opt.value();

	std::cout << "wav sample rate: " << wave.samples_per_sec << "\n";

	uint64_t sample_count =
	    ((wave.samples.size() - 1) / wave.samples_per_sec) * 48000;

	/*	std::vector<float> one_channel_sample;

	  for (int i = 0; i < wave.samples.size(); i += wave.channels) {
	    one_channel_sample.push_back(wave.samples[i]);
	  }

	  std::vector<float> resampleSampels = resampleAudio(
	      one_channel_sample, wave.samples_per_sec, wave.samples_per_sec *
	  4, 1);

	  uint64_t increment = (resampleSampels.size() - 1) * 0.1;

	  std::cout << "samplerate new: " << wave.samples_per_sec * 4 << "\n";

	  for (int i = 0; i < resampleSampels.size() - 1; i += increment + 1) {
	    std::ofstream csv(std::format("output{}.csv", i));
	    csv << "idx;value\n";
	    std::cout << "samples " << i << " to " << i + increment << "\n";
	    std::vector<float> samples;
	    for (int j = i; j < i + increment; j++) {
	      samples.push_back(resampleSampels[j]);
	    }

	    std::vector<std::complex<double>> fft_results = fft(samples);
	    std::vector<double> results;

	    double max_result = 0;

	    for (std::complex<double> result_complex : fft_results) {
	      double result = std::abs(result_complex);
	      if (result > max_result) {
	        max_result = result;
	      }
	      results.push_back(result);
	    }

	    for (int j = 0; j < results.size(); j++) {
	      results[j] = results[j] / max_result;
	    }

	    for (int j = 0; j < (results.size() - 1) / 2; j++) {
	      double value = results[j];
	      if (value < 0.5)
	        continue;
	      csv << j << ";" << value << "\n";
	    }
	    csv.close();
	  }*/
	float orig_freq = 51.91;
	// float new_freq = 25.96;
	// float new_freq = 103.83;
	float new_freq = 25.96;

	float one_Hz_samplerate = wave.samples_per_sec / new_freq;

	std::vector<float> pitch_ajusted =
	    resampleAudio(wave.samples, wave.samples_per_sec,
	                  one_Hz_samplerate * orig_freq, wave.channels);
	// return 0;

	Pa_Initialize();
	PaStream *stream;
	PaStreamParameters outputParameters;

	int device = -1;

	// std::cout << "argument count: " << argc << "\n";
	for (int i = 1; i < argc; i++) {
		// std::cout << i << ": " << argv[i] << "\n";
		if (strcmp(argv[i], "-d ")) {
			if (argc > i + 1) {
				i++;
				device = atoi(argv[i]);
				std::cout << "device: " << device << "\n";
			}
		}
	}

	if (device == -1) {
		for (int i = 0; i < Pa_GetDeviceCount(); i++) {
			const PaDeviceInfo *pInfo = Pa_GetDeviceInfo(i);
			if (pInfo->maxOutputChannels < 1) {
				std::cout << i << ": output device doesn't have channels"
				          << "\n";
				continue;
			}
			std::cout << i << ": " << pInfo->name
			          << ", hostApi:" << pInfo->hostApi << "\n";
		}
		std::cout << "pick audio device: ";
		std::cin >> device;
	}

	outputParameters.device = device;
	if (outputParameters.device == paNoDevice) {
		return -1;
	}

	const PaDeviceInfo *pInfo = Pa_GetDeviceInfo(device);
	if (pInfo != 0) {
		printf("Output device name: '%s' samplerate: %f\n", pInfo->name,
		       pInfo->defaultSampleRate);
		sample_rate = pInfo->defaultSampleRate;
	} else {
		sample_rate = sample_rate_default;
	}

	std::cout << "sample rate: " << sample_rate << "\n";

	samples = resampleAudio(pitch_ajusted, wave.samples_per_sec, sample_rate,
	                        wave.channels);

	sample_offset = 1.0 /
	                (wave.channels * (ceilf(wave.bits_per_sample / 2.0))) *
	                (wave.samples_per_sec / sample_rate);
	float seconds =
	    (samples.size() / (float)wave.channels) / (float)sample_rate;

	outputParameters.channelCount = 2; /* stereo output */
	outputParameters.sampleFormat =
	    paFloat32; /* 32 bit floating point output */
	outputParameters.suggestedLatency =
	    Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
	outputParameters.hostApiSpecificStreamInfo = NULL;

	PaError err = Pa_OpenStream(&stream, NULL, /* no input */
	                            &outputParameters, sample_rate, 0,
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
		std::cout << "playing seconds: " << seconds << "\n";
		Pa_Sleep(seconds * 1000);
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

	if (K <= 1) {
		return input;
	}

	std::vector<std::complex<double>> even(N / 2);
	std::vector<std::complex<double>> odd(N / 2);

	for (int i = 0; 2 * i < N; i++) {
		even.push_back(input[2 * i]);
		odd.push_back(input[2 * i + 1]);
	}

	std::vector<std::complex<double>> even_output = fft(even);
	std::vector<std::complex<double>> odd_output = fft(odd);
	std::vector<std::complex<double>> output(N);
	for (int k = 0; k < K / 2; k++) {
		std::complex<double> t =
		    std::polar(1.0, 2 * M_PI * k / N) * odd_output[k];
		output[k] = even_output[k] + t;
		output[N / 2 + k] = even_output[k] - t;
	}

	return output;
}
