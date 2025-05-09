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
#include <sstream>
#include <sys/types.h>
#include <variant>
#include <vector>

#define PRINT_MIDI_INFO 0
#define envelope_log 0

uint32_t progress = 0;
std::vector<float> samples_to_play;
bool looped = false;
int loop_start;
int loop_end;
std::ofstream info_log_file("log.txt");
float volume = 0.75;
bool isFinished = false;

int paCallback(const void *inputBuffer, void *outputBuffer,
               unsigned long framesPerBuffer,
               const PaStreamCallbackTimeInfo *timeInfo,
               PaStreamCallbackFlags statusFlags, void *userData) {
	float *out = (float *)outputBuffer;
	for (int i = 0; i < framesPerBuffer; i++) {
		if (!isFinished) {
			*out++ = samples_to_play[progress] * volume;
			*out++ = samples_to_play[progress + 1] * volume;
			progress += 2;
			if (progress >= samples_to_play.size()) {
				isFinished = true;
				std::cout << "complete\n";
			}
		} else {
			*out++ = 0;
			*out++ = 0;
		}
	}
	if (isFinished) {
		return paComplete;
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
	float pan = 0;
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
	for (uint32_t i = 0; i < this->samples.size(); i++) {
		if (this->samples[i].low_key <= note &&
		    note <= this->samples[i].high_key) {
			sample = samples[i];
			found = true;
			break;
		}
	}

	if (!found) {
		std::cout << "sample of note " << (uint32_t)note << " not found\n";
		return frames;
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

	float volume_left = (0.5 - sample.pan) * volume;
	float volume_right = (0.5 + sample.pan) * volume;

	uint32_t index = 0;
	uint32_t size = 0;
	if (sample.samplemodes == 0 || sample.samplemodes == 2) {
		for (uint32_t i = 0; i < frames_needed; i++) {
			float sample_value = 0;
			if (index < sample.samples.size()) {
				sample_value = sample.samples[index];
			}

			frames.push_back(sample_value * volume_left);  // left
			frames.push_back(sample_value * volume_right); // right

			index += 1;
		}
	} else if (sample.samplemodes == 1 || sample.samplemodes == 3) {
		for (uint32_t i = 0; i < frames_needed; i++) {
			float sample_value = sample.samples[index];
			frames.push_back(sample_value * volume_left);  // left
			frames.push_back(sample_value * volume_right); // right

			index += 1;
			if (index >= sample.loop_end - 1) {
				float before = sample.samples[index];
				index = sample.loop_start;
				float after = sample.samples[index];
			}
		}
	}

	frames = resampleAudio(frames, sample.sample_rate, new_sample_rate, 2);

	frames = resampleAudio(frames, sample.sample_rate, sample_rate, 2);

	float volume_envelope = sample.sustainVol;
	for (uint32_t i = 0; i < frames.size(); i += 2) {
		float time_ms = ((i / 2) / sample_rate) * 1000.0f;
		volume_envelope = sample.sustainVol;
		if (time_ms < sample.attackVolTime) {
			volume_envelope = time_ms / sample.attackVolTime;
#if envelope_log
			std::cout << time_ms << " ";
			std::cout << sample.attackVolTime
			          << " attack| vol: " << volume_envelope * 100 << "%\n";
#endif
		} else if (time_ms < sample.attackVolTime + sample.holdVolTime) {
			volume_envelope = 1.0f;
#if envelope_log
			std::cout << time_ms << " ";
			std::cout << sample.attackVolTime + sample.holdVolTime
			          << " hold| vol: " << volume_envelope * 100 << "%\n";
#endif
		} else if (time_ms < sample.attackVolTime + sample.holdVolTime +
		                         sample.decayVolTime) {
			float percentage =
			    (time_ms - sample.attackVolTime - sample.holdVolTime) /
			    sample.decayVolTime;
			volume_envelope = 1 - ((1.0f - sample.sustainVol) * percentage);
#if envelope_log
			std::cout << time_ms << " ";
			std::cout << sample.attackVolTime + sample.holdVolTime +
			                 sample.decayVolTime
			          << " decay| percent: " << volume_envelope * 100 << "%\n";
#endif
		}
		frames[i] *= volume_envelope;     // left
		frames[i + 1] *= volume_envelope; // right
	}

	std::vector<float> release_frames;

	if (sample.samplemodes == 1) {
		for (uint32_t i = 0; i < frames_needed_release; i++) {
			float sample_value = 0;
			if (index < sample.samples.size()) {
				sample_value = sample.samples[index];
			}

			release_frames.push_back(sample_value * volume_left);  // left
			release_frames.push_back(sample_value * volume_right); // right

			index += 1;
			if (index >= sample.loop_end - 1) {
				float before = sample.samples[index];
				index = sample.loop_start;
				float after = sample.samples[index];
				// frames.push_back(before + (after - before) * 0.5);
			}
		}
	} else if (sample.samplemodes == 3) {
		for (uint32_t i = 0; i < frames_needed_release; i++) {
			float sample_value = 0;
			if (index < sample.samples.size()) {
				sample_value = sample.samples[index];
			}

			release_frames.push_back(sample_value * volume_left);  // left
			release_frames.push_back(sample_value * volume_right); // right

			index += 1;
		}
	} else {
		for (uint32_t i = 0; i < frames_needed_release; i++) {
			float sample_value = 0;
			if (index < sample.samples.size()) {
				sample_value = sample.samples[index];
			}

			release_frames.push_back(sample_value * volume_left);  // left
			release_frames.push_back(sample_value * volume_right); // right

			index += 1;
		}
	}

	release_frames =
	    resampleAudio(release_frames, sample.sample_rate, new_sample_rate, 2);

	release_frames =
	    resampleAudio(release_frames, sample.sample_rate, sample_rate, 2);

	frames.reserve(release_frames.size());
	for (uint32_t i = 0; i < release_frames.size(); i += 2) {
		float percent = (float)(i / 2.0f) / (release_frames.size() / 2.0f);
		float volume_release = volume_envelope - (volume_envelope * percent);
#if envelope_log
		std::cout << "release percent: " << percent * 100
		          << "% vol: " << volume_release * 100 << "%\n";
#endif
		frames.push_back(release_frames[i] * volume_release);     // left
		frames.push_back(release_frames[i + 1] * volume_release); // right
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
	if (frames.size() % 2 != 0) {
		std::cout << "frame count not even\n";
		exit(99);
	}

	return frames;
}

std::vector<float> Preset::get_audio(uint8_t note, float seconds,
                                     float sample_rate, float volume) {
	std::vector<float> raw_samples =
	    instrument.get_audio(note, seconds, sample_rate, volume);
	std::vector<float> delayBuffer_left;
	std::vector<float> delayBuffer_right;
	std::vector<float> samples;
	float decayFactor = reverbEffect / 100.0f;
	delayBuffer_left.resize(500, 0);
	delayBuffer_right.resize(500, 0);
	samples.reserve(raw_samples.size());
	if (decayFactor == 0.0f) {
		samples = raw_samples;
	} else {
		for (uint32_t i = 0; i < raw_samples.size(); i += 2) {
			float sample_left = raw_samples[i];
			float sample_right = raw_samples[i + 1];
			delayBuffer_left.push_back(sample_left);
			delayBuffer_right.push_back(sample_right);

			for (uint32_t j = 0; j < delayBuffer_left.size(); j++) {
				delayBuffer_left[j] *= decayFactor;
				delayBuffer_right[j] *= decayFactor;
			}
			float output_left = 0.0;
			float output_right = 0.0;
			for (uint32_t j = 0; j < delayBuffer_left.size(); j++) {
				output_left += delayBuffer_left[j];
				output_right += delayBuffer_right[j];
			}

			delayBuffer_left.erase(delayBuffer_left.begin());
			delayBuffer_right.erase(delayBuffer_right.begin());

			samples.push_back(output_left);
			samples.push_back(output_right);
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
		Sample sample;
		uint32_t key = -1;
		for (uint16_t i = instrument_bag_min.inst_gen_ndx;
		     i < instrument_bag_max_ex.inst_gen_ndx && !end; i++) {
			sf_inst_gen generator = sf.instrument_generator_list[i];
			switch (generator.oper) {
			case sf_generator::keyRange: {
				sf_ranges_type keys = generator.amount.ranges;
				info_log_file << "key range: " << (uint32_t)keys.low << "-"
				              << (uint32_t)keys.high << "\n";
				sample.low_key = keys.low;
				sample.high_key = keys.high;
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
				sample.pan = generator.amount.sh_amount / 1000;
				break;
			}
			case sf_generator::sampleModes: {
				info_log_file << "samplemodes: " << generator.amount.w_amount
				              << "\n";
				sample.samplemodes = generator.amount.w_amount;
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
				sample.attackVolTime =
				    pow(2.0f, (generator.amount.sh_amount / 1200.0f)) * 1000;

				break;
			}
			case sf_generator::holdVolEnv: {
				info_log_file
				    << "holdVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms|orig: " << generator.amount.sh_amount << "\n";
				sample.holdVolTime =
				    pow(2.0f, (generator.amount.sh_amount / 1200.0f)) * 1000;
				break;
			}
			case sf_generator::decayVolEnv: {
				info_log_file
				    << "decayVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms|orig: " << generator.amount.sh_amount << "\n";
				sample.decayVolTime =
				    pow(2.0f, (generator.amount.sh_amount / 1200.0f)) * 1000;
				break;
			}
			case sf_generator::sustainVolEnv: {
				info_log_file
				    << "sustainVolEnv: "
				    << std::pow(10.0f, (-generator.amount.sh_amount / 200.0f))
				    << "f|orig: " << generator.amount.sh_amount << "\n";
				sample.sustainVol = std::pow(
				    10.0f, ((-generator.amount.sh_amount / 10.0f) / 20.0f));
				break;
			}
			case sf_generator::releaseVolEnv: {
				info_log_file
				    << "releaseVolEnv: "
				    << pow(2.0f, (generator.amount.sh_amount / 1200.0f)) *
				           1000.0f
				    << " ms|orig: " << generator.amount.sh_amount << "\n";
				sample.releaseVolTime =
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
		sample.loop_start = sfsample.start_loop - sfsample.start;
		sample.loop_end = sfsample.end_loop - sfsample.start;
		sample.sample_rate = sfsample.sample_rate;
		if (key == -1) {
			sample.original_key = sfsample.original_key;
		} else {
			sample.original_key = key;
		}
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
			bank[sf.preset_headers[i].preset] = preset;
		}
	}

	return bank;
}

struct Note {
	uint32_t start;
	uint32_t duration;
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
	uint32_t start;
	uint8_t amount;
};

struct Track {
	std::string name;
	int preset;
	std::vector<std::variant<Note, Controller>> events;
};

void process_midi_event(Track &track, uint32_t absolute_time,
                        float ms_per_midiclock, MIDI_MIDI_EVENT midi) {
	uint8_t type = (0xF0 & midi.type) >> 4;
	uint8_t channel = 0x0F & midi.type;
	if (channel != 0) {
		std::cout
		    << "this program expects only midi channel \"0\" to be used\n";
		exit(-1);
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
			if (note.midi_note == key && note.duration == 0) {
				note.duration = absolute_time - note.start;
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
			Note note = {.start = absolute_time,
			             .duration = 0,
			             .midi_note = key,
			             .velocity = velocity};
			track.events.push_back(note);
		} else {
			for (int32_t i = track.events.size() - 1; i >= 0; i--) {
				if (!std::holds_alternative<Note>(track.events[i])) {
					continue;
				}
				Note note = std::get<Note>(track.events[i]);
				if (note.midi_note == key && note.duration == 0) {
					note.duration = absolute_time - note.start;
#if PRINT_MIDI_INFO
					std::cout << "note: " << (uint32_t)note.midi_note
					          << " duration: " << note.duration_ticks << "\n";
#endif
					track.events[i] = note;
					return;
				}
			}
#if PRINT_MIDI_INFO
			std::cout << "note on not found\n";
#endif
		}
	} break;
	case 0xB: {
		Controller controller = {.type = (Controller_type)midi.data[0],
		                         .start = absolute_time,
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
		std::cerr << "type: 0x" << std::hex << std::uppercase
		          << (uint32_t)midi.type << std::dec << std::nouppercase
		          << "\n";
	}
}

void write_envelopes_to_csvs(std::array<Preset, 128> bank) {
	for (uint32_t i = 0; i < bank.size(); i++) {
		Preset preset = bank[i];
		for (uint32_t j = 0; j < preset.instrument.samples.size(); j++) {
			Sample sample_source = preset.instrument.samples[j];
			uint32_t smm = 1000;
			std::vector<float> samples;
			samples.resize(smm * 50, 1.0f);
			uint32_t loop_end = samples.size() - 10;
			Sample sample = {
			    .samples = samples,
			    .loop_start = 10,
			    .loop_end = loop_end,
			    .sample_rate = smm,
			    .original_key = 60,
			    .low_key = 0,
			    .high_key = 127,
			    .attackVolTime = sample_source.attackVolTime,
			    .holdVolTime = sample_source.holdVolTime,
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

			std::vector<float> env_samples =
			    inst.get_audio(sample_source.low_key, 10, smm, 1.0f);

			std::stringstream string_stream_left;
			string_stream_left << "./output/envelope/waveform_envelope";
			string_stream_left << i << ".." << j << "L";
			string_stream_left << ".csv";
			std::cout << string_stream_left.str() << "\n";

			std::ofstream waveform1(string_stream_left.str());
			waveform1 << "x,y\n";
			for (uint32_t i = 0; i < env_samples.size(); i += 2) {
				float sample = env_samples[i];
				waveform1 << i / 2 << "," << sample << "\n";
			}
			waveform1.close();

			std::stringstream string_stream_right;
			string_stream_right << "./output/envelope/waveform_envelope";
			string_stream_right << i << ".." << j << "R";
			string_stream_right << ".csv";
			std::cout << string_stream_right.str() << "\n";

			std::ofstream waveform2(string_stream_right.str());
			waveform2 << "x,y\n";
			for (uint32_t i = 1; i < env_samples.size(); i += 2) {
				float sample = env_samples[i];
				waveform2 << floor(i / 2.0f) << "," << sample << "\n";
			}
			waveform2.close();
		}
	}
}

int main(int argc, char **argv) {
	float sample_rate;
	float sample_rate_default = 48000.0;

	Pa_Initialize();
	int device = -1;
	int track_index = -1;
	std::string sf_file_name = "MUS_LAST_BOSS.sf2";
	std::string mid_file_name = "MUS_LAST_BOSS.mid";

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-c") == 0) {
			if (argc > i + 1) {
				i++;
				device = atoi(argv[i]);
				std::cout << "device: " << device << "\n";
			}
		} else if (strcmp(argv[i], "-d") == 0) {
			device = Pa_GetDefaultOutputDevice();
			std::cout << "device: " << device << "\n";
		} else if (strcmp(argv[i], "-t") == 0) {
			if (argc > i + 1) {
				i++;
				track_index = atoi(argv[i]);
				std::cout << "track_index: " << track_index << "\n";
			}
		} else if (strcmp(argv[i], "-M") == 0) {
			if (argc > i + 1) {
				i++;
				mid_file_name = argv[i];
				std::cout << "midi file: " << mid_file_name << "\n";
			}
		} else if (strcmp(argv[i], "-S") == 0) {
			if (argc > i + 1) {
				i++;
				sf_file_name = argv[i];
				std::cout << "sf2 file: " << sf_file_name << "\n";
			}
		}
	}

	SF2File sf2_file(sf_file_name);
	sf_sample sample_header = sf2_file.sample_headers[40];
	std::array<Preset, 128> bank = get_bank(sf2_file, 0);
	/*
	    Preset preset = bank[70];

	    return 0;
	*/
	MIDI_FILE midi_file;
	std::optional<MIDI_FILE> midi_file_opt =
	    read_midi_file(mid_file_name.c_str());

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
#if PRINT_MIDI_INFO
				std::cout << "SYSEX\n";
				std::cout << "data: 0x" << std::hex;
#endif
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
					track.name = std::string(data.begin(), data.end());
					std::cout << "text: " << data.data() << "\n";
				}
				if (meta.type == 0x51) {
					uint32_t time = meta.data[0];
					time = time << 8 | meta.data[1];
					time = time << 8 | meta.data[2];
#if PRINT_MIDI_INFO
					std::cout << "time: " << time << "\n";
#endif
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

	PaStream *stream;
	PaStreamParameters outputParameters;

	float seconds = 40;

	if (device == -1) {
		for (int i = 0; i < Pa_GetDeviceCount(); i++) {
			const PaDeviceInfo *pInfo = Pa_GetDeviceInfo(i);
			if (pInfo->maxOutputChannels < 1) {
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

	float mult = (us_per_midiclock /
	              (midi_file.header.division.ticks_per_quarter_note)) /
	             1000000;

	int track_min = 0;
	int track_max = tracks.size();
	if (track_index != -1) {
		track_min = track_index;
		track_max = track_index + 1;
	}
#if 1
	int32_t events_count = 0;
	for (uint32_t i = track_min; i < track_max; i++) {
		for (uint32_t j = 0; j < tracks[i].events.size(); j++) {
			events_count += 1;
		}
	}
	uint32_t events_done = 0;
	// 3 6
	int allowed[] = {0, 1, 2, 4, 5, 7, 8, 9, 10, 11, 12};
	for (uint32_t i = track_min; i < track_max; i++) {
		/*if (std::find(std::begin(allowed), std::end(allowed), i) ==
		    std::end(allowed)) {
		    continue;
		}*/
		Preset preset = bank[tracks[i].preset];
		Instrument instrument = preset.instrument;
		float expression = 1.0f;
		float pan = 0.0f;
		std::vector<float> track_samples;
		for (uint32_t j = 0; j < tracks[i].events.size(); j++) {
			std::cout << events_done + 1 << "/" << events_count << "\n";
			if (std::holds_alternative<Note>(tracks[i].events[j])) {
				Note note = std::get<Note>(tracks[i].events[j]);
				float volume = expression * (note.velocity / 127.0f);
				float start_sec = note.start * mult;
				int start_index = std::round(sample_rate * start_sec) * 2;
				float duration_sec = note.duration * mult;
				// std::cout << "sec: " << start_sec << " dur: " << duration_sec
				//          << "\n";
				std::vector<float> samples = preset.get_audio(
				    note.midi_note, duration_sec, sample_rate, volume);
				uint32_t max_index = samples.size() - 1 + start_index;
				if (max_index >= track_samples.size()) {
					track_samples.resize(max_index + 1);
				}
				for (uint32_t index = 0; index < samples.size(); index += 2) {
					if (pan == 0) {
						track_samples[start_index + index] +=
						    samples[index] * volume; // left
						track_samples[start_index + index + 1] +=
						    samples[index + 1] * volume; // right
					} else {
						track_samples[start_index + index] +=
						    samples[index] * volume * (0.5 - pan); // left
						track_samples[start_index + index + 1] +=
						    samples[index + 1] * volume * (0.5 + pan); // right
					}
				}
			} else if (std::holds_alternative<Controller>(
			               tracks[i].events[j])) {
				Controller controller =
				    std::get<Controller>(tracks[i].events[j]);
				if (controller.type == Controller_type::Expression) {
					float prev_expression = expression;
					expression = controller.amount / 127.0f;
					float start_sec = controller.start * mult;
					int start_index = std::round(sample_rate * start_sec) * 2;
					if (start_index < track_samples.size()) {
						for (uint32_t index = start_index;
						     index < track_samples.size(); index++) {

							track_samples[index] *=
							    expression / prev_expression;
						}
					}
					// std::cout << "val: " << (uint32_t)controller.amount
					//           << ", express: " << std::fixed
					//           << std::setprecision(10)
					//           << (uint32_t)(expression * 100) << "\n";
				} else if (controller.type == Controller_type::Pan) {
					float prev_pan = pan;
					if (controller.amount != 64) { // not center
						pan = (controller.amount / 127.0) * 2 - 1;
					} else { // center
						pan = 0;
					}
					float start_sec = controller.start * mult;
					int start_index = std::round(sample_rate * start_sec) * 2;
					if (pan != 0) {
						if (start_index < track_samples.size()) {
							for (uint32_t index = start_index;
							     index + 1 < track_samples.size(); index += 2) {
								track_samples[index] *=
								    (0.5 - pan) / (0.5 - prev_pan);
								track_samples[index + 1] *=
								    (0.5 + pan) / (0.5 + prev_pan);
							}
						}
					} else {
						if (start_index < track_samples.size()) {
							for (uint32_t index = start_index;
							     index + 1 < track_samples.size(); index += 2) {
								track_samples[index] /= (0.5 - prev_pan);
								track_samples[index + 1] /= (0.5 + prev_pan);
							}
						}
					}
					// std::cout << "pan: " << (uint32_t)controller.amount
					//           << "%\n";
				}
			}
			events_done += 1;
		}
		float max_sample =
		    *std::max_element(track_samples.begin(), track_samples.end());
		if (max_sample > 1.0f) {
			for (auto &sample : track_samples) {
				sample /= max_sample;
			}
		}
		if (track_samples.size() > samples_to_play.size()) {
			samples_to_play.resize(track_samples.size());
		}
		for (uint32_t j = 0; j < track_samples.size(); j++) {
			samples_to_play[j] += track_samples[j];
		}
	}
	uint32_t zeros = 0;
	bool nonzero_found = false;
	if (samples_to_play.size() == 0) {
		std::cout << "empty\n";
		exit(1000);
	}

	for (uint32_t i = samples_to_play.size() - 1; i >= 0 && !nonzero_found;
	     i--) {
		if (samples_to_play[i] == 0) {
			zeros += 1;
		} else {
			nonzero_found = true;
		}
	}
	samples_to_play.erase(samples_to_play.end() - zeros, samples_to_play.end());

	/*zeros = 0;
	nonzero_found = false;
	std::cout << "beforeloop " << zeros << " " << !nonzero_found << "\n";
	for (uint32_t i = 0; i < samples_to_play.size() && !nonzero_found; i++) {
	    if (samples_to_play[i] == 0) {
	        zeros += 1;
	    } else {
	        nonzero_found = true;
	    }
	}
	samples_to_play.erase(samples_to_play.begin(),
	                      samples_to_play.begin() + zeros);*/
#else
	write_envelopes_to_csvs(bank);
	return 0;
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

	uint32_t smm = 48000;
	;
	std::vector<float> samples;
	samples.resize(smm * 50, 1.0f);
	Sample sample = {
	    .samples = samples,
	    .sample_rate = 22025, //(uint32_t)floor(smm * 1.4),
	    .original_key = 60,
	    .low_key = 0,
	    .high_key = 127,
	    .attackVolTime = sample_source.attackVolTime,
	    .holdVolTime = sample_source.holdVolTime,
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

	outputParameters.channelCount = 2; /* stereo output */
	outputParameters.sampleFormat =
	    paFloat32; /* 32 bit floating point output */
	outputParameters.suggestedLatency =
	    Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
	outputParameters.hostApiSpecificStreamInfo = NULL;
#if 0
	std::ofstream waveform2("./output/waveform_audio.csv");
	std::ofstream waveform3("./output/waveform_audio.bin", std::ios::binary);

	waveform2 << "index,sample\n";
	for (uint32_t i = 0; i < samples_to_play.size(); i += 2) {
		waveform2 << i / 2 << "," << samples_to_play[i] << "\n";
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
	std::cout << "length: " << (float)samples_to_play.size() / 2 / sample_rate
	          << "s\n";
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
		for (uint32_t i = 0;
		     i < ((samples_to_play.size() / 2.0f) / sample_rate) * 10; i++) {
			float seconds = (i + 1) / 10.0;
			Pa_Sleep(100);
			std::cout << std::fixed << std::setprecision(1)
			          << "current: " << seconds << "/"
			          << ((samples_to_play.size() / 2.0f) / sample_rate) << "\n"
			          << std::defaultfloat;
		}
		std::cout << "done\n";
		Pa_StopStream(stream);
	}
	Pa_CloseStream(stream);
	Pa_Terminate();

	return 0;
}
