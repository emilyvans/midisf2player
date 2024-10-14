/** @file paex_sine.c
    @ingroup examples_src
    @brief Play a sine wave for several seconds.
    @author Ross Bencina <rossb@audiomulch.com>
    @author Phil Burk <philburk@softsynth.com>
*/
/*
 * $Id: paex_sine.c 1752 2011-09-08 03:21:55Z philburk $
 *
 * This program uses the PortAudio Portable Audio Library.
 * For more information see: http://www.portaudio.com/
 * Copyright (c) 1999-2000 Ross Bencina and Phil Burk
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files
 * (the "Software"), to deal in the Software without restriction,
 * including without limitation the rights to use, copy, modify, merge,
 * publish, distribute, sublicense, and/or sell copies of the Software,
 * and to permit persons to whom the Software is furnished to do so,
 * subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
 * ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * The text above constitutes the entire PortAudio license; however,
 * the PortAudio community also makes the following non-binding requests:
 *
 * Any person wishing to distribute modifications to the Software is
 * requested to send the modifications to the original developer so that
 * they can be incorporated into the canonical version. It is also
 * requested that these non-binding requests be included along with the
 * license above.
 */
#include "portaudio.h"
#include <cmath>
#include <fcntl.h>
#include <math.h>
#include <pulse/timeval.h>
#include <stdio.h>
#include <unistd.h>

#define NUM_SECONDS (15)
// #define SAMPLE_RATE (44100)
#define SAMPLE_RATE (48000)
#define FRAMES_PER_BUFFER (256)

#ifndef M_PI
#define M_PI (3.14159265)
#endif

#define FREQ 220

#define TABLE_SIZE (FREQ)

class Sine {
public:
	Sine() : stream(0), left_phase(0), right_phase(0) {
		FILE *f = fopen("sin.csv", "w");
		/* initialise sinusoidal wavetable */
		double phaseIncrement = (FREQ * M_PI * 2.0) / SAMPLE_RATE;
		for (int i = 0; i < TABLE_SIZE; i++) {
			sine[i] = (float)sin(2.0 * M_PI * ((double)i) / (double)TABLE_SIZE);
			fprintf(f, "%d,%f\n", i + 1, sine[i]);
		}
		fclose(f);
		sprintf(message, "No Message");
		printf("%lu", paFramesPerBufferUnspecified);
	}

	bool open(PaDeviceIndex index) {
		PaStreamParameters outputParameters;

		outputParameters.device = index;
		if (outputParameters.device == paNoDevice) {
			return false;
		}

		const PaDeviceInfo *pInfo = Pa_GetDeviceInfo(index);
		if (pInfo != 0) {
			printf("Output device name: '%s'\r", pInfo->name);
		}

		outputParameters.channelCount = 2; /* stereo output */
		outputParameters.sampleFormat =
				paFloat32; /* 32 bit floating point output */
		outputParameters.suggestedLatency =
				Pa_GetDeviceInfo(outputParameters.device)->defaultLowOutputLatency;
		outputParameters.hostApiSpecificStreamInfo = NULL;

		PaError err = Pa_OpenStream(
				&stream, NULL, /* no input */
				&outputParameters, SAMPLE_RATE, FRAMES_PER_BUFFER,
				paClipOff, /* we won't output out of range samples so don't bother
		                    clipping them */
				&Sine::paCallback, this /* Using 'this' for userData so we can cast
		                                 to Sine* in paCallback method */
		);

		if (err != paNoError) {
			/* Failed to open stream to device !!! */
			return false;
		}

		err = Pa_SetStreamFinishedCallback(stream, &Sine::paStreamFinished);

		if (err != paNoError) {
			Pa_CloseStream(stream);
			stream = 0;

			return false;
		}

		return true;
	}

	bool close() {
		if (stream == 0)
			return false;

		PaError err = Pa_CloseStream(stream);
		stream = 0;

		return (err == paNoError);
	}

	bool start() {
		if (stream == 0)
			return false;
		PaError err = Pa_StartStream(stream);

		return (err == paNoError);
	}

	bool stop() {
		if (stream == 0)
			return false;

		PaError err = Pa_StopStream(stream);

		return (err == paNoError);
	}

private:
	/* The instance callback, where we have access to every method/variable in
	 * object of class Sine */
	int paCallbackMethod(const void *inputBuffer, void *outputBuffer,
	                     unsigned long framesPerBuffer,
	                     const PaStreamCallbackTimeInfo *timeInfo,
	                     PaStreamCallbackFlags statusFlags) {
		float *out = (float *)outputBuffer;
		unsigned long i;

		(void)timeInfo; /* Prevent unused variable warnings. */
		(void)statusFlags;
		(void)inputBuffer;

		for (i = 0; i < framesPerBuffer; i++) {
			double deltaFrame = (double)prevSteps / (double)SAMPLE_RATE;
			float f1 = SinWave(130.813, deltaFrame);
			float f2 = SinWave(164.814, deltaFrame);
			float f3 = SinWave(195.998, deltaFrame);
			float combination = f1 + f2 + f3;
			*out++ = combination; /* left */
			*out++ = combination; /* right */
			left_phase += 1;
			if (left_phase >= TABLE_SIZE)
				left_phase -= TABLE_SIZE;
			right_phase += 2; /* higher pitch so we can distinguish left and right. */
			if (right_phase >= TABLE_SIZE)
				right_phase -= TABLE_SIZE;
			prevSteps++;
		}

		return paContinue;
	}

	float SinWave(float freq, double delta) {
		return sin(freq * 2.0 * M_PI * delta);
	}

	/* This routine will be called by the PortAudio engine when audio is needed.
	** It may called at interrupt level on some machines so don't do anything
	** that could mess up the system like calling malloc() or free().
	*/
	static int paCallback(const void *inputBuffer, void *outputBuffer,
	                      unsigned long framesPerBuffer,
	                      const PaStreamCallbackTimeInfo *timeInfo,
	                      PaStreamCallbackFlags statusFlags, void *userData) {
		/* Here we cast userData to Sine* type so we can call the instance
		   method paCallbackMethod, we can do that since we called Pa_OpenStream
		   with 'this' for userData */
		return ((Sine *)userData)
		    ->paCallbackMethod(inputBuffer, outputBuffer, framesPerBuffer, timeInfo,
		                       statusFlags);
	}

	void paStreamFinishedMethod() { printf("Stream Completed: %s\n", message); }

	/*
	 * This routine is called by portaudio when playback is done.
	 */
	static void paStreamFinished(void *userData) {
		return ((Sine *)userData)->paStreamFinishedMethod();
	}

	PaStream *stream;
	int prevSteps = 0;
	float sine[TABLE_SIZE];
	int left_phase;
	int right_phase;
	char message[20];
};

class ScopedPaHandler {
public:
	ScopedPaHandler() : _result(Pa_Initialize()) {}
	~ScopedPaHandler() {
		if (_result == paNoError) {
			Pa_Terminate();
		}
	}

	PaError result() const { return _result; }

private:
	PaError _result;
};

/*******************************************************************/
int main(void);
int main(void) {
	Sine sine;

	printf("PortAudio Test: output sine wave. SR = %d\n", SAMPLE_RATE);

	ScopedPaHandler paInit;
	if (paInit.result() != paNoError)
		goto error;

	if (sine.open(Pa_GetDefaultOutputDevice())) {
		if (sine.start()) {
			printf("Play for %d seconds.\n", NUM_SECONDS);
			Pa_Sleep(NUM_SECONDS * 1000);

			sine.stop();
		}

		sine.close();
	}

	printf("Test finished.\n");
	return paNoError;

error:
	fprintf(stderr, "An error occurred while using the portaudio stream\n");
	fprintf(stderr, "Error number: %d\n", paInit.result());
	fprintf(stderr, "Error message: %s\n", Pa_GetErrorText(paInit.result()));
	return 1;
}
