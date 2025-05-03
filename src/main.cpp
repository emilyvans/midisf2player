#include "midi_file.hpp"
#include "midi_note_freq.hpp"
#include "sf2_file.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <ostream>
#include <portaudio.h>
#include <ratio>
#include <sys/types.h>
#include <variant>
#include <vector>

uint32_t progress = 0;
std::vector<float> samples_to_play;
bool looped = false;
int loop_start;
int loop_end;
std::ofstream info_log_file("log.txt");
float volume = 0.75;

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
			return paAbort;
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
	float attackVolTime = 0.0f;
	float holdVolTime = 0.0f;
	float decayVolTime = 0.0f;
	float sustainVol = 1.0f;
	float releaseVolTime = 0.0f;
	uint32_t samplemodes = 0;
};

class Instrument {
	friend Instrument get_instrument(SF2File, int);
	friend struct Preset;

  private:
  public:
	std::string name;
	std::vector<Sample> samples;
	std::vector<float> get_audio(uint8_t note, float seconds, float sample_rate,
	                             float volume);
};

struct Preset {
	Instrument instrument;
	float reverbEffect = 0.0f;
	std::vector<float> get_audio(uint8_t note, float seconds, float sample_rate,
	                             float volume);
};

std::vector<float> Instrument::get_audio(uint8_t note, float seconds,
                                         float sample_rate, float volume) {
	std::vector<float> frames;

	Sample sample;
	bool found = false;
	std::cout << name << "\n";
	for (uint32_t i = 0; i < this->samples.size(); i++) {
		// std::cout << (samples[i].low_key <= note) << "-"
		//           << (note <= samples[i].high_key) << "\n";
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
	uint32_t frames_needed_release = sample.sample_rate *
	                                 (sample.releaseVolTime / 1000) *
	                                 (newfreq / origfreq);
	frames.reserve(frames_needed);

	uint32_t index = 0;
	uint32_t size = 0;
	std::cout << (uint32_t)sample.samplemodes << "\n";
	if (sample.samplemodes == 0 || sample.samplemodes == 2) {
		std::cout << "noloop\n";
		for (uint32_t i = 0; i < frames_needed; i++) {
			// frames.push_back(sample.samples[index] * volume);

			float sample_value = 0;
			if (index < sample.samples.size()) {
				sample_value = sample.samples[index];
			}

			float time_ms = (i / sample_rate) * 1000.0f;

			frames.push_back(sample_value * volume);

			index += 1;
		}

		size = frames.size();
	} else if (sample.samplemodes == 1 || sample.samplemodes == 3) {
		std::cout << "loop\n";
		for (uint32_t i = 0; i < frames_needed; i++) {
			// frames.push_back(sample.samples[index] * volume);

			float sample_value = sample.samples[index];
			frames.push_back(sample_value * volume);

			index += 1;
			if (index >= sample.loop_end - 1) {
				float before = sample.samples[index];
				index = sample.loop_start;
				float after = sample.samples[index];
				// frames.push_back(before + (after - before) * 0.5);
			}
		}
		size = frames.size();
	}

	frames = resampleAudio(frames, sample.sample_rate, new_sample_rate, 1);

	frames = resampleAudio(frames, sample.sample_rate, sample_rate, 1);

	float volume_envelope = sample.sustainVol;
	for (uint32_t i = 0; i < frames.size(); i++) {
		float time_ms = (i / sample_rate) * 1000.0f;
		volume_envelope = sample.sustainVol;
		if (time_ms < sample.attackVolTime) {
			volume_envelope = time_ms / sample.attackVolTime;
			// std::cout << time_ms << " ";
			// std::cout << sample.attackVolTime
			//           << " attack| vol: " << volume_envelope * 100 << "%\n";
		} else if (time_ms < sample.attackVolTime + sample.holdVolTime) {
			volume_envelope = 1.0f;
			// std::cout << time_ms << " ";
			// std::cout << sample.attackVolTime + sample.holdVolTime
			//           << " hold| vol: " << volume_envelope * 100 << "%\n";
		} else if (time_ms < sample.attackVolTime + sample.holdVolTime +
		                         sample.decayVolTime) {
			float percentage =
			    (time_ms - sample.attackVolTime - sample.holdVolTime) /
			    sample.decayVolTime;
			// volume_envelope = (sample.sustainVol) * (1 - percentage);
			volume_envelope = 1 - ((1.0f - sample.sustainVol) * percentage);
			// std::cout << time_ms << " ";
			/// std::cout << sample.attackVolTime + sample.holdVolTime +
			//                 sample.decayVolTime
			//          << " decay| percent: " << volume_envelope * 100 <<
			//          "%\n";
		}
		frames[i] *= volume_envelope;
	}

	std::vector<float> release_frames;

	if (sample.samplemodes == 1) {
		for (uint32_t i = 0; i < frames_needed_release; i++) {
			float sample_value = 0;
			float volume_release =
			    volume_envelope -
			    (volume_envelope * i / (frames_needed_release - 1));
			if (index < sample.samples.size()) {
				sample_value = sample.samples[index];
			}

			release_frames.push_back(sample_value * volume_release * volume);

			index += 1;
			if (index >= sample.loop_end - 1) {
				float before = sample.samples[index];
				index = sample.loop_start;
				float after = sample.samples[index];
				// frames.push_back(before + (after - before) * 0.5);
			}
		}
	} else if (sample.samplemodes == 3) {
		std::cout << "release\n";
		// index = sample.loop_end;
		for (uint32_t i = 0; i < frames_needed_release; i++) {
			float sample_value = 0;
			float volume_release =
			    volume_envelope -
			    (volume_envelope * i / (frames_needed_release - 1));
			if (index < sample.samples.size()) {
				sample_value = sample.samples[index];
			}

			release_frames.push_back(sample_value * volume_release * volume);

			index += 1;
		}
	} else {
		for (uint32_t i = 0; i < frames_needed_release; i++) {
			float sample_value = 0;
			float volume_release =
			    volume_envelope -
			    (volume_envelope * i / (frames_needed_release - 1));
			if (index < sample.samples.size()) {
				sample_value = sample.samples[index];
			}

			release_frames.push_back(sample_value * volume_release * volume);

			index += 1;
		}
	}

	release_frames =
	    resampleAudio(release_frames, sample.sample_rate, new_sample_rate, 1);

	release_frames =
	    resampleAudio(release_frames, sample.sample_rate, sample_rate, 1);

	frames.reserve(release_frames.size());
	for (uint32_t i = 0; i < release_frames.size(); i++) {
		frames.push_back(release_frames[i]);
	};

	/*if (!((frames.size() / sample_rate) <
	          (seconds + sample.releaseVolTime / 1000) + 0.01 &&
	      (frames.size() / sample_rate) >
	          (seconds + sample.releaseVolTime / 1000) - 0.01)) {
	    std::cout << "too long: " << (frames.size() / sample_rate)
	              << "!=" << (seconds + sample.releaseVolTime / 1000)
	              << " -+0.01\n";
	    exit(100);
	}*/

	return frames;
}

std::vector<float> Preset::get_audio(uint8_t note, float seconds,
                                     float sample_rate, float volume) {
	std::vector<float> raw_samples =
	    instrument.get_audio(note, seconds, sample_rate, volume);
	std::vector<float> delayBuffer;
	std::vector<float> samples;
	float decayFactor = reverbEffect / 100.0f;
	delayBuffer.resize(500, 0);
	samples.reserve(raw_samples.size());
	if (decayFactor == 0.0f) {
		std::cout << "zero: " << decayFactor << "\n";
		samples = raw_samples;
	} else {
		for (uint32_t i = 0; i < raw_samples.size(); i++) {
			// std::cout << "ds: " << i << "/" << raw_samples.size() << "\n";
			float sample = raw_samples[i];
			delayBuffer.push_back(sample);

			for (uint32_t j = 0; j < delayBuffer.size(); j++) {
				// std::cout << "db1: " << j << "/" << delayBuffer.size() <<
				// "\n";
				delayBuffer[j] *= decayFactor;
			}
			float output = 0.0;
			for (uint32_t j = 0; j < delayBuffer.size(); j++) {
				// std::cout << "db2: " << j << "/" << delayBuffer.size() <<
				// "\n";
				output += delayBuffer[j];
			}

			delayBuffer.erase(delayBuffer.begin());

			samples.push_back(output);
		}
	}
	float max = 0.0f;

	for (float sample : samples) {
		if (std::abs(sample) > max) {
			max = std::abs(sample);
		}
	}

	for (uint32_t i = 0; i < samples.size(); i++) {
		samples[i] /= max;
	}

	return samples;
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
		uint32_t samplemodes = 0;
		float attackVolTime = 0.0f;
		float holdVolTime = 0.0f;
		float decayVolTime = 0.0f;
		float sustainVol = 1.0f;
		float releaseVolTime = 0.0f;
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
				samplemodes = generator.amount.w_amount;
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
				    << " ms|orig: " << generator.amount.sh_amount << "\n";
				attackVolTime =
				    pow(2.0f, (generator.amount.sh_amount / 1200.0f)) * 1000;

				break;
			}
			case sf_generator::holdVolEnv: {
				info_log_file
				    << "holdVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms|orig: " << generator.amount.sh_amount << "\n";
				holdVolTime =
				    pow(2.0f, (generator.amount.sh_amount / 1200.0f)) * 1000;
				break;
			}
			case sf_generator::decayVolEnv: {
				info_log_file
				    << "decayVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms|orig: " << generator.amount.sh_amount << "\n";
				decayVolTime =
				    pow(2.0f, (generator.amount.sh_amount / 1200.0f)) * 1000;
				break;
			}
			case sf_generator::sustainVolEnv: {
				info_log_file
				    << "sustainVolEnv: "
				    << std::pow(10.0f, (-generator.amount.sh_amount / 200.0f))
				    << "f|orig: " << generator.amount.sh_amount << "\n";
				sustainVol = std::pow(
				    10.0f, ((-generator.amount.sh_amount / 10.0f) / 20.0f));
				break;
			}
			case sf_generator::releaseVolEnv: {
				info_log_file
				    << "releaseVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms|orig: " << generator.amount.sh_amount << "\n";
				releaseVolTime =
				    pow(2.0f, (generator.amount.sh_amount / 1200.0f)) * 1000;
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
				exit(1);
				break;
			}
		}
		sf_sample sfsample = sf.sample_headers[sample_number];
		info_log_file << "sample name: " << sfsample.sample_name << "\n";
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
		sample.attackVolTime = attackVolTime;
		sample.holdVolTime = holdVolTime;
		sample.decayVolTime = decayVolTime;
		sample.sustainVol = sustainVol;
		sample.releaseVolTime = releaseVolTime;
		sample.samplemodes = samplemodes;
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
			float percentage = 0.1f * generator.amount.sh_amount;
			preset.reverbEffect = percentage;
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
			exit(999);
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

std::array<Preset, 128> get_bank(SF2File sf, uint32_t bank_number) {
	std::array<Preset, 128> bank = {};

	// skips last element because its a teminal record which is complety zeros
	for (uint32_t i = 0; i < sf.preset_headers.size() - 1; i++) {
		if (sf.preset_headers[i].bank == bank_number) {
			Preset preset = getPreset(sf, i);
			// std::cout << i + 1 << "/" << sf.preset_headers.size() - 1 <<
			// "\n";
			bank[sf.preset_headers[i].preset] = preset;
			// bank.push_back(preset);
		}
	}
	info_log_file << 54
	              << "|loop start:" << bank[54].instrument.samples[0].loop_start
	              << "\n";
	info_log_file << 54
	              << "|loop end:" << bank[54].instrument.samples[0].loop_end
	              << "\n";

	return bank;
}

struct Note {
	uint32_t start_ticks;
	//	float start_ms;
	uint32_t duration_ticks;
	//	float duration_ms;
	uint8_t midi_note;
	uint32_t velocity;
};

enum class Controller_type : uint8_t {
	Balance = 0x08,
	Pan = 0x0A,
	Expression = 0x0B,
};

struct Controller {
	Controller_type type;
	uint8_t amount;
};

struct Track {
	int preset;
	std::vector<std::variant<Note, Controller>> events;
};

#define PRINT_MIDI_INFO 1

void process_midi_event(Track &track, uint32_t absolute_time,
                        float ms_per_midiclock, MIDI_MIDI_EVENT midi) {
	uint8_t type = (0xF0 & midi.type) >> 4;
	uint8_t channel = 0x0F & midi.type;
	if (channel != 0) {
		exit(9999);
	}
	switch (type) {
	case 0x8: {
		uint8_t key = midi.data[0];
		uint8_t velocity = midi.data[1];
#if PRINT_MIDI_INFO
		std::cout << "off note: " << (uint32_t)key
		          << " vel: " << (uint32_t)velocity
		          << " timestamp: " << absolute_time << "\n";
#endif
		for (int32_t i = track.events.size() - 1; i >= 0; i--) {
			if (!std::holds_alternative<Note>(track.events[i])) {
				continue;
			}
			Note note = std::get<Note>(track.events[i]);
			if (note.midi_note == key && note.duration_ticks == 0) {
				note.duration_ticks = absolute_time - note.start_ticks;
#if PRINT_MIDI_INFO
				std::cout << "note: " << (uint32_t)note.midi_note
				          << " duration: " << note.duration_ticks << "\n";
#endif
				track.events[i] = note;
				return;
			}
		}
		std::cout << "note on not found\n";
	} break;
	case 0x9: {
		uint8_t key = midi.data[0];
		uint8_t velocity = midi.data[1];
#if PRINT_MIDI_INFO
		std::cout << "on note: " << (uint32_t)key
		          << " vel: " << (uint32_t)velocity
		          << " timestamp: " << absolute_time << "\n";
#endif
		if (velocity != 0) {
			// std::cout << "noteon\n";
			Note note = {.start_ticks = absolute_time,
			             .duration_ticks = 0,
			             .midi_note = key,
			             .velocity = velocity};
			// std::cout << track.notes.size() << "\n";
			track.events.push_back(note);
			// std::cout << track.notes.size() << "\n";
		} else {
			for (int32_t i = track.events.size() - 1; i >= 0; i--) {
				if (!std::holds_alternative<Note>(track.events[i])) {
					continue;
				}
				Note note = std::get<Note>(track.events[i]);
				if (note.midi_note == key && note.duration_ticks == 0) {
					note.duration_ticks = absolute_time - note.start_ticks;
#if PRINT_MIDI_INFO
					std::cout << "note: " << (uint32_t)note.midi_note
					          << " duration: " << note.duration_ticks << "\n";
#endif
					track.events[i] = note;
					return;
				}
			}
			std::cout << "note on not found\n";
		}
	} break;
	case 0xB: {
		Controller controller = {.type = (Controller_type)midi.data[0],
		                         .amount = midi.data[1]};
#if PRINT_MIDI_INFO
		std::cout << "control change:\n";
		std::cout << "type: " << (uint32_t)controller.type
		          << ", amount: " << (uint32_t)controller.amount << "\n";
#endif
		track.events.push_back(controller);
	} break;
	case 0xC: {
		track.preset = midi.data[0];
#if PRINT_MIDI_INFO
		std::cout << "program: " << (uint32_t)midi.data[0] << "\n";
#endif
	} break;
	default:
		std::cout << "type: 0x" << std::hex << std::uppercase
		          << (uint32_t)midi.type << std::dec << std::nouppercase
		          << "\n";
		exit(999);
	}
}

int main(int argc, char **argv) {
	float sample_rate;
	float sample_rate_default = 48000.0;

	SF2File sf2_file("./MUS_LAST_BOSS.sf2");
	sf_sample sample_header = sf2_file.sample_headers[40];
	std::array bank = get_bank(sf2_file, 0);

	MIDI_FILE midi_file;
	std::optional<MIDI_FILE> midi_file_opt =
	    read_midi_file("./MUS_LAST_BOSS.mid");
	// read_midi_file("./test2.mid");

	if (!midi_file_opt.has_value()) {
		std::cout << "couldn't read midi file: \"" << "./MUS_LAST_BOSS.mid"
		          << "\"\n";
		exit(1);
	}

	midi_file = midi_file_opt.value();
	float us_per_midiclock = 0;
	std::vector<Track> tracks;
	float length = 0;

	for (int i = 0; i < midi_file.tracks.size(); i++) {
		MIDI_TRACK midi_track = midi_file.tracks[i];
		Track track;
		uint32_t absolute_time = 0;
		for (int j = 0; j < midi_track.events.size(); j++) {
			MIDI_EVENT event = midi_track.events[j];
			absolute_time += event.delta_time;
			if (std::holds_alternative<MIDI_SYSEX_EVENT>(event.event)) {
				MIDI_SYSEX_EVENT sysex =
				    std::get<MIDI_SYSEX_EVENT>(event.event);
				std::cout << "SYSEX\n";
				std::cout << "data: 0x" << std::hex;
				for (int index = 0; index < sysex.data.size(); index++) {
					std::cout << (uint32_t)sysex.data[index];
				}
				std::cout << std::dec << "\n";
			} else if (std::holds_alternative<MIDI_MIDI_EVENT>(event.event)) {
				MIDI_MIDI_EVENT midi = std::get<MIDI_MIDI_EVENT>(event.event);
				// std::cout << "MIDI\n";
				process_midi_event(track, absolute_time, us_per_midiclock,
				                   midi);
			} else if (std::holds_alternative<MIDI_META_EVENT>(event.event)) {
				MIDI_META_EVENT meta = std::get<MIDI_META_EVENT>(event.event);
				if (meta.type == 0x2F) {
					continue;
				}
				// std::cout << "META\n";
				// std::cout << "meta type: 0x" << std::hex <<
				// (uint32_t)meta.type
				//           << std::dec << "\n";
				if (meta.type == 0x21) {
					// std::cout << "port: " << (uint32_t)meta.data[0] << "\n";
				}
				if (meta.type == 0x3) {
					std::vector<uint8_t> data = meta.data;
					data.push_back(0);
					// std::cout << "text: " << data.data() << "\n";
				}
				if (meta.type == 0x51) {
					uint32_t time = meta.data[0];
					time = time << 8 | meta.data[1];
					time = time << 8 | meta.data[2];
					std::cout << "time: " << time << "\n";
					us_per_midiclock = (float)time;
				}
			}
		}
		std::cout << "--------------------" << i + 1 << "-----------------\n";
		if (length < absolute_time) {
			length = absolute_time;
		}
		tracks.push_back(track);
	}

	// return 0;
	Pa_Initialize();
	PaStream *stream;
	PaStreamParameters outputParameters;

	float seconds = 40;

	int device = -1;
	int track_index = -1;

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
		} else if (strcmp(argv[i], "-t") == 0) {
			if (argc > i + 1) {
				i++;
				track_index = atoi(argv[i]);
				std::cout << "track_index: " << track_index << "\n";
			}
		}
	}

	if (device == -1) {
		for (int i = 0; i < Pa_GetDeviceCount(); i++) {
			const PaDeviceInfo *pInfo = Pa_GetDeviceInfo(i);
			if (pInfo->maxOutputChannels < 1) {
				// std::cout << i << ": output device doesn't have channels"
				//           << "\n";
				continue;
			}
			const PaHostApiInfo *hInfo = Pa_GetHostApiInfo(pInfo->hostApi);
			std::cout << i << ": " << pInfo->name << ", hostApi:" << hInfo->name
			          << "\n";
		}
		std::cout << "pick audio device: ";
		std::cin >> device;
		std::cin.ignore();
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

	// std::cout << bank[9].instrument.name << "\n";

	// samples_to_play = bank[11].instrument.get_audio(21, seconds,
	// sample_rate);
	//
	float mult = //(1 / length) * 93;
	    (us_per_midiclock /
	     (midi_file.header.division.ticks_per_quarter_note)) /
	    1000000;

	int track_min = 0;
	int track_max = tracks.size();
	if (track_index != -1) {
		track_min = track_index;
		track_max = track_index + 1;
	}
#if 1
	for (uint32_t i = track_min; i < track_max; i++) {
		Preset preset = bank[tracks[i].preset];
		Instrument instrument = preset.instrument;
		float expression = 1.0f;
		// std::reverse(tracks[i].notes.begin(), tracks[i].notes.end());
		for (uint32_t j = 0; j < tracks[i].events.size(); j++) {
			std::cout << i + 1 << "/" << tracks.size() << "| " << j + 1 << "/"
			          << tracks[i].events.size() << "\n";
			if (std::holds_alternative<Note>(tracks[i].events[j])) {
				Note note = std::get<Note>(tracks[i].events[j]);
				float volume = expression * (note.velocity / 127.0f);
				float start_sec = note.start_ticks * mult;
				int start_index = std::round(sample_rate * start_sec);
				float duration_sec = note.duration_ticks * mult;
				std::cout << "sec: " << start_sec << " dur: " << duration_sec
				          << "\n";
				std::vector<float> samples = preset.get_audio(
				    note.midi_note, duration_sec, sample_rate, volume);
				uint32_t max_index = samples.size() - 1 + start_index;
				if (max_index >= samples_to_play.size()) {
					samples_to_play.resize(max_index + 1);
				}
				std::cout << samples.size() << "\n";
				for (uint32_t index = 0; index < samples.size(); index++) {
					float original_sample =
					    samples_to_play[start_index + index];
					float new_sample = samples[index] * volume;

					samples_to_play[start_index + index] += new_sample;
				}
			} else if (std::holds_alternative<Controller>(
			               tracks[i].events[j])) {
				Controller controller =
				    std::get<Controller>(tracks[i].events[j]);
				if (controller.type == Controller_type::Expression) {
					expression = controller.amount / 127.0f;
					std::cout << "val: " << (uint32_t)controller.amount
					          << ", express: " << std::fixed
					          << std::setprecision(10)
					          << (uint32_t)(expression * 100) << "\n";
				}
			}
		}
	}
#else
	Preset preset = bank[track_index];
	int note = 60;
	int sec = 20;
	assert(sec < 50);
	samples_to_play = preset.get_audio(note, sec, sample_rate, 1.0f);
	std::cout << "song size: " << samples_to_play.size() << "\n";
	if (samples_to_play.size() != 0) {
		std::cout << "middle sample: "
		          << samples_to_play[(samples_to_play.size() - 1) / 2] << "\n";
	}

	Sample sample_source;
	bool found_source;
	for (Sample sample : preset.instrument.samples) {
		std::cout << (uint32_t)sample.low_key << " <= " << note
		          << " <= " << (uint32_t)sample.high_key << "\n";
		if (sample.low_key <= note && note <= sample.high_key) {
			sample_source = sample;
			found_source = true;
			break;
		}
	}
	if (!found_source) {
		std::cout << "not found\n";
		exit(9999);
	}

	uint32_t smm = 1000;
	std::vector<float> samples;
	samples.resize(smm * 50, 1.0f);
	Sample sample = {
	    .samples = samples,
	    .sample_rate = smm,
	    .original_key = 60,
	    .low_key = 0,
	    .high_key = 127,
	    .attackVolTime = sample_source.attackVolTime, // 1000,
	    .holdVolTime = sample_source.holdVolTime,     // 0.00000602385,
	    .decayVolTime = sample_source.decayVolTime,
	    .sustainVol = sample_source.sustainVol,
	    .releaseVolTime = sample_source.releaseVolTime,
	    .samplemodes = sample_source.samplemodes,
	};
	std::vector<Sample> sm;
	sm.push_back(sample);
	Instrument inst = {
	    .name = "t",
	    .samples = sm,
	};

	std::vector<float> env_samples = inst.get_audio(note, sec, smm, 1.0f);

	std::ofstream waveform1("./output/waveform_envelope.csv");
	waveform1 << "x,y\n";
	for (uint32_t i = 0; i < env_samples.size(); i++) {
		float sample = env_samples[i];
		waveform1 << i << "," << sample << "\n";
	}
	waveform1.close();
#endif
	// return 0;

	outputParameters.channelCount = 2; /* stereo output */
	outputParameters.sampleFormat =
	    paFloat32; /* 32 bit floating point output */
	outputParameters.suggestedLatency =
	    Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
	outputParameters.hostApiSpecificStreamInfo = NULL;
#if 1
	std::ofstream waveform2("./output/waveform_audio.csv");
	std::ofstream waveform3("./output/waveform_audio.bin", std::ios::binary);

	waveform2 << "index,sample\n";
	for (uint32_t i = 0; i < samples_to_play.size(); i++) {
		waveform2 << i << "," << samples_to_play[i] << "\n";
		waveform3 << samples_to_play[i];
	}
	waveform2.close();
	waveform3.close();

#endif

	float max_sample =
	    *std::max_element(samples_to_play.begin(), samples_to_play.end());
	if (max_sample > 1.0f) {
		for (auto &sample : samples_to_play) {
			sample /= max_sample;
		}
	}
	std::cout << length << "\n"
	          << midi_file.header.division.ticks_per_quarter_note << "\n";
	// return 0;
	std::cout << "start\n";
	std::cin.get();

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
		for (uint32_t i = 0; i < (samples_to_play.size() / sample_rate) * 10;
		     i++) {
			float seconds = (i + 1) / 10.0;
			std::cout << std::fixed << std::setprecision(1)
			          << "current: " << seconds << "/"
			          << (samples_to_play.size() / sample_rate) << "\n"
			          << std::defaultfloat;
			Pa_Sleep(100);
		}
		// Pa_Sleep(length * 1000);
		std::cout << "done\n";
		Pa_StopStream(stream);
	}
	Pa_CloseStream(stream);
	Pa_Terminate();

	return 0;
}
