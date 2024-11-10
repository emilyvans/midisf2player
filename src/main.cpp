#include "sf2_file.hpp"
#include <cstdint>
#include <cstring>
#include <iostream>
#include <portaudio.h>
#include <vector>

uint32_t progress = 0;
std::vector<float> samples;
bool looped = false;
int loop_start;
int loop_end;

int paCallback(const void *inputBuffer, void *outputBuffer,
               unsigned long framesPerBuffer,
               const PaStreamCallbackTimeInfo *timeInfo,
               PaStreamCallbackFlags statusFlags, void *userData) {
	float *out = (float *)outputBuffer;
	for (int i = 0; i < framesPerBuffer; i++) {
		float sample;
		if ((!looped && progress == samples.size()) ||
		    (looped && progress == loop_end)) {
			looped = true;
			int pre_progress = progress;
			progress = loop_start + 1;
			sample = 0.5 * samples[loop_start] + 0.5 * samples[pre_progress];

		} else {
			sample = samples[progress];
		}
		*out++ = sample * 0.5;
		*out++ = sample * 0.5;
		progress++;
	}

	return paContinue;
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
	float sample_rate;
	float sample_rate_default = 48000.0;
	std::cout << "Hello, World!\n";
	SF2File sf2_file("./MUS_LAST_BOSS.sf2");
	sf_sample sample_header = sf2_file.sample_headers[9];

	std::vector<float> pre_samples;

	for (int i = sample_header.start; i < sample_header.end + 1; i++) {
		pre_samples.push_back(sf2_file.samples[i]);
	}

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

	loop_start = (sample_header.start_loop - sample_header.start) / 22050.0 *
	             sample_rate;
	loop_end =
	    (sample_header.end_loop - sample_header.start) / 22050.0 * sample_rate;

	std::cout << "sample rate: " << sample_rate << "\n";
	std::cout << "start-loop: " << loop_start << ", end-loop: " << loop_end
	          << "\n";

	samples = resampleAudio(pre_samples, 22050, sample_rate, 1);

	float seconds = 20;

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
	                            paCallback, nullptr);

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

	return 0;
}
