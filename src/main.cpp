#include "midi_note_freq.hpp"
#include "sf2_file.hpp"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <ostream>
#include <portaudio.h>
#include <vector>

uint32_t progress = 0;
std::vector<float> samples_to_play;
bool looped = false;
int loop_start;
int loop_end;
std::ofstream info_log_file("log.txt");
float volume = 0.25;

int paCallback(const void *inputBuffer, void *outputBuffer,
               unsigned long framesPerBuffer,
               const PaStreamCallbackTimeInfo *timeInfo,
               PaStreamCallbackFlags statusFlags, void *userData) {
	float *out = (float *)outputBuffer;
	for (int i = 0; i < framesPerBuffer; i++) {
		float sample;
#if 0
		if ((!looped && progress == samples.size()) ||
		    (looped && progress == loop_end)) {
			looped = true;
			int pre_progress = progress;
			progress = loop_start + 1;
			// sample = 0.5 * samples[loop_start] + 0.5 * samples[pre_progress];
			sample = samples[progress];
		} else {
			sample = samples[progress];
		}
#else
		sample = samples_to_play[progress];
		// std::cout << "running: " << progress << "/" << samples_to_play.size()
		//           << "\n";
#endif
		*out++ = sample * volume;
		*out++ = sample * volume;
		progress++;
		if (progress >= samples_to_play.size()) {
			std::cout << "complete\n";
			return paComplete;
		}
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

struct Sample {
	std::vector<float> samples;
	uint32_t loop_start;
	uint32_t loop_end;
	uint32_t sample_rate;
	uint8_t original_key;
	uint8_t low_key;
	uint8_t high_key;
};

class Instrument {
	friend Instrument get_instrument(SF2File, int);

  private:
  public:
	std::string name;
	std::vector<Sample> samples;
	std::vector<float> get_audio(uint8_t note, float seconds,
	                             float sample_rate);
};

struct Preset {
	Instrument instrument;
};

std::vector<float> Instrument::get_audio(uint8_t note, float seconds,
                                         float sample_rate) {
	std::vector<float> frames;

	Sample sample;
	bool found = false;
	std::cout << this->samples.size() << "\n";
	for (uint32_t i = 0; i < this->samples.size(); i++) {
		std::cout << (samples[i].low_key <= note) << "-"
		          << (note <= samples[i].high_key) << "\n";
		if (this->samples[i].low_key <= note &&
		    note <= this->samples[i].high_key) {
			sample = samples[i];
			found = true;
			break;
		}
	}

	if (!found) {
		std::cout << "sample of note " << (uint32_t)note << " not found\n";
		exit(1);
	}

	double origfreq = midi_note_freq[sample.original_key];
	double newfreq = midi_note_freq[note];
	float new_sample_rate = sample.sample_rate / newfreq * origfreq;

	uint32_t frames_needed =
	    sample.sample_rate * seconds * (newfreq / origfreq);
	frames.reserve(frames_needed);

	uint32_t index = 0;
	for (uint32_t i = 0; i < frames_needed; i++) {
		frames.push_back(sample.samples[index]);
		index += 1;
		if (index >= sample.loop_end) {
			index = sample.loop_start;
		}
	}

	frames = resampleAudio(frames, sample.sample_rate, new_sample_rate, 1);

	frames = resampleAudio(frames, sample.sample_rate, sample_rate, 1);

	return frames;
}

Instrument get_instrument(SF2File sf, int index) {
	sf_inst sf_instrument = sf.instrument_names_indices[index];
	Instrument instrument;
	uint32_t instrument_zone_amount =
	    sf.instrument_names_indices[index + 1].inst_bag_ndx -
	    sf_instrument.inst_bag_ndx;

	instrument.name = sf_instrument.inst_name;

	info_log_file << "instrument        |" << sf_instrument.inst_name << "\n";
	info_log_file << "instrument bag ndx|" << sf_instrument.inst_bag_ndx
	              << "\n";
	info_log_file << "instrument bag cnd|" << instrument_zone_amount << "\n";

	for (int i = 0; i < instrument_zone_amount; i++) {
		sf_inst_bag instrument_bag_min =
		    sf.instument_index_list[sf_instrument.inst_bag_ndx + i];
		sf_inst_bag instrument_bag_max_ex =
		    sf.instument_index_list[sf_instrument.inst_bag_ndx + i + 1];
		info_log_file << "instrument mod idx|"
		              << instrument_bag_min.inst_mod_ndx << " "
		              << instrument_bag_max_ex.inst_mod_ndx -
		                     instrument_bag_min.inst_mod_ndx
		              << " " << instrument_bag_min.inst_mod_ndx << "-"
		              << instrument_bag_max_ex.inst_mod_ndx << "\n";
		info_log_file << "instrument gen idx|"
		              << instrument_bag_min.inst_gen_ndx << " "
		              << instrument_bag_max_ex.inst_gen_ndx -
		                     instrument_bag_min.inst_gen_ndx
		              << " " << instrument_bag_min.inst_gen_ndx << "-"
		              << instrument_bag_max_ex.inst_gen_ndx << "\n";

		bool end = false;
		uint32_t sample_number;
		uint32_t low_key = 0;
		uint32_t high_key = 127;
		uint32_t key = -1;
		for (uint16_t i = instrument_bag_min.inst_gen_ndx;
		     i < instrument_bag_max_ex.inst_gen_ndx && !end; i++) {
			sf_inst_gen generator = sf.instrument_generator_list[i];
			switch (generator.oper) {
			case sf_generator::keyRange: {
				sf_ranges_type keys = generator.amount.ranges;
				info_log_file << "key range: " << (uint32_t)keys.low << "-"
				              << (uint32_t)keys.high << "\n";
				low_key = keys.low;
				high_key = keys.high;
				break;
			}
			case sf_generator::velRange: {
				sf_ranges_type vels = generator.amount.ranges;
				info_log_file << "velocity range: " << (uint32_t)vels.low << "-"
				              << (uint32_t)vels.high << "\n";
				break;
			}
			case sf_generator::initialAttenuation: {
				info_log_file << "initial Attenuation: "
				              << generator.amount.sh_amount / 10 << "dB\n";
				break;
			}
			case sf_generator::pan: {
				info_log_file << "pan: " << generator.amount.sh_amount / 10
				              << "%\n";
				break;
			}
			case sf_generator::sampleModes: {
				info_log_file << "samplemodes: " << generator.amount.w_amount
				              << "\n";
				break;
			}
			case sf_generator::overridingRootKey: {
				info_log_file
				    << "overridingRootKey: " << generator.amount.sh_amount
				    << "\n";
				key = generator.amount.sh_amount;
				break;
			}
			case sf_generator::attackVolEnv: {
				info_log_file
				    << "attackVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms\n";
				break;
			}
			case sf_generator::holdVolEnv: {
				info_log_file
				    << "holdVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms\n";
				break;
			}
			case sf_generator::decayVolEnv: {
				info_log_file
				    << "decayVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms\n";
				break;
			}
			case sf_generator::sustainVolEnv: {
				info_log_file
				    << "sustainVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms\n";
				break;
			}
			case sf_generator::releaseVolEnv: {
				info_log_file
				    << "releaseVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms\n";
				break;
			}
			case sf_generator::sampleID: {
				info_log_file << "sampleID: " << generator.amount.w_amount
				              << "\n";
				sample_number = generator.amount.w_amount;
				end = true;
				break;
			}
			default:
				std::cout << "generator: \"" << generator.oper
				          << "\" is not supported\n";
				break;
			}
		}
		sf_sample sfsample = sf.sample_headers[sample_number];
		Sample sample;
		sample.loop_start = sfsample.start_loop - sfsample.start;
		sample.loop_end = sfsample.end_loop - sfsample.start;
		sample.sample_rate = sfsample.sample_rate;
		if (key == -1) {
			sample.original_key = sfsample.original_key;
		} else {
			sample.original_key = key;
		}
		sample.high_key = high_key;
		sample.low_key = low_key;
		for (int sample_index = sfsample.start; sample_index < sfsample.end;
		     sample_index++) {
			sample.samples.push_back(sf.samples[sample_index]);
		}
		instrument.samples.push_back(sample);
	}
	return instrument;
}

Preset getPreset(SF2File sf, uint32_t index) {
	Preset preset;
	sf_preset_header header = sf.preset_headers[index];
	uint32_t bag_amount =
	    sf.preset_headers[index + 1].perset_bag_ndx - header.perset_bag_ndx;

	if (header.preset == 0) {
		std::cout << "couldn't find preset of index: " << index << "\n";
	}

	info_log_file << "preset            |" << header.preset_name << "\n";
	info_log_file << "preset nr         |" << header.preset << "\n";
	info_log_file << "preset bank       |" << header.bank << "\n";
	info_log_file << "preset bag index  |" << header.perset_bag_ndx << "\n";
	info_log_file << "preset library    |" << header.library << "\n";
	info_log_file << "preset genre      |" << header.genre << "\n";
	info_log_file << "preset morphology |" << header.morphology << "\n";
	info_log_file << "bag amount        |" << bag_amount << "\n";
	if (bag_amount != 1) {
		std::cout << "the preset has more or less then one bag. this program "
		             "does not "
		             "support more than one\n";
		exit(1);
	}
	sf_preset_bag preset_indices_min =
	    sf.preset_index_list[header.perset_bag_ndx];
	sf_preset_bag preset_indices_max_ex =
	    sf.preset_index_list[header.perset_bag_ndx + 1];
	info_log_file << "modulators        |"
	              << preset_indices_max_ex.ModNdx - preset_indices_min.ModNdx
	              << " " << preset_indices_min.ModNdx << "-"
	              << preset_indices_max_ex.ModNdx << "\n";

	info_log_file << "generators        |"
	              << preset_indices_max_ex.GenNdx - preset_indices_min.GenNdx
	              << " " << preset_indices_min.GenNdx << "-"
	              << preset_indices_max_ex.GenNdx << "\n";

	if ((preset_indices_max_ex.ModNdx - preset_indices_min.ModNdx) != 0) {
		std::cout << "this program doesn't support modulators\n";
	}

	bool end = false;
	int32_t instrument_number = -1;
	for (uint16_t i = preset_indices_min.GenNdx;
	     i < preset_indices_max_ex.GenNdx && !end; i++) {
		sf_preset_gen generator = sf.preset_generator_list[i];
		switch (generator.oper) {
		case sf_generator::reverbEffectsSend: {
			float percentage = 0.1 * generator.amount.sh_amount;
			if (percentage != 0) {
				info_log_file << "reverbEffect: " << percentage << "%\n";
			}
			break;
		}
		case sf_generator::instrument: {
			info_log_file << "instument: " << generator.amount.w_amount << "\n";
			instrument_number = generator.amount.w_amount;
			end = true;
			break;
		}
		default:
			std::cout << "generator: \"" << generator.oper
			          << "\" is not supported\n";
			break;
		}
	}
	// instruments
	Instrument instrument = get_instrument(sf, instrument_number);
	for (int i = 0; i < instrument.samples.size(); i++) {
		info_log_file << i << "|loop start:" << instrument.samples[i].loop_start
		              << "\n";
		info_log_file << i << "|loop end:" << instrument.samples[i].loop_end
		              << "\n";
		info_log_file << i
		              << "|sample rate:" << instrument.samples[i].sample_rate
		              << "\n";
		info_log_file << i << "|original key:"
		              << (uint32_t)instrument.samples[i].original_key << "\n";
		info_log_file << i << "|high key:"
		              << (uint32_t)instrument.samples[i].high_key << "\n";
		info_log_file << i
		              << "|low key:" << (uint32_t)instrument.samples[i].low_key
		              << "\n";
	}

	preset.instrument = instrument;

	return preset;
}

std::vector<Preset> get_bank(SF2File sf, uint32_t bank_number) {
	std::vector<Preset> bank;

	// skips last element because its a teminal record which is complety zeros
	for (uint32_t i = 0; i < sf.preset_headers.size() - 1; i++) {
		if (sf.preset_headers[i].bank == bank_number) {
			Preset preset = getPreset(sf, i);
			bank.push_back(preset);
		}
	}
	info_log_file << 9
	              << "|loop start:" << bank[9].instrument.samples[0].loop_start
	              << "\n";
	info_log_file << 9 << "|loop end:" << bank[9].instrument.samples[0].loop_end
	              << "\n";

	return bank;
}

int main(int argc, char **argv) {
	float sample_rate;
	float sample_rate_default = 48000.0;
	SF2File sf2_file("./MUS_LAST_BOSS.sf2");
	sf_sample sample_header = sf2_file.sample_headers[40];
	std::vector<Preset> bank = get_bank(sf2_file, 0);
	// 1, 2, 30, 38, 45, 47, 48, 52, 53, 54, 55, 56(fav), 81, 126
	/*
	    for (int i = 0; i < sf2_file.preset_headers.size(); i++) {
	        sf_preset_header preset_header = sf2_file.preset_headers[i];
	        std::cout << "--------------------------------\n";
	        // preset header
	        std::cout << "preset:            " << preset_header.preset_name <<
	   "\n"; std::cout << "preset_nr:         " << preset_header.preset << "\n";
	        std::cout << "preset_bank:       " << preset_header.bank << "\n";
	        std::cout << "preset_bag_index:  " << preset_header.perset_bag_ndx
	                  << "\n";
	        std::cout << "preset_library:    " << preset_header.library << "\n";
	        std::cout << "preset_genre:      " << preset_header.genre << "\n";
	        std::cout << "preset_morphology: " << preset_header.morphology <<
	   "\n";

	        // preset generator
	        // std::cout << "generator: " << i << "\n";
	        // std::cout << "operator:  " << (int)preset_gen.oper << "\n";
	        // std::cout << "amount:    " << preset_gen.amount.w_amount << "\n";
	    }
	    std::cout << "--------------------------------\n";
	*/

	Pa_Initialize();
	PaStream *stream;
	PaStreamParameters outputParameters;

	int device = -1;

	// std::cout << "argument count: " << argc << "\n";
	for (int i = 1; i < argc; i++) {
		// std::cout << i << ": |" << argv[i] << "|\n";
		if (strcmp(argv[i], "-c") == 0) {
			if (argc > i + 1) {
				i++;
				device = atoi(argv[i]);
				std::cout << "device: " << device << "\n";
			}
		} else if (strcmp(argv[i], "-d") == 0) {
			device = Pa_GetDefaultOutputDevice();
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
			const PaHostApiInfo *hInfo = Pa_GetHostApiInfo(pInfo->hostApi);
			std::cout << i << ": " << pInfo->name << ", hostApi:" << hInfo->name
			          << "\n";
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

	float seconds = 40;

	std::cout << bank[9].instrument.name << "\n";

	samples_to_play = bank[11].instrument.get_audio(21, seconds, sample_rate);

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
