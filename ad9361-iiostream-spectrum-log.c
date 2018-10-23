/*
 * libiio - AD9361 IIO streaming example
 *
 * Copyright (C) 2014 IABG mbH
 * Author: Michael Feilen <feilen_at_iabg.de>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 **/

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <signal.h>
#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <unistd.h>
#include <errno.h>
#include <pthread.h>

#ifdef __APPLE__
#include <iio/iio.h>
#else
#include <iio.h>
#endif

/* helper macros */
#define MHZ(x) ((long long)(x*1000000.0 + .5))
#define GHZ(x) ((long long)(x*1000000000.0 + .5))

#define AVERAGES 1;
#define FREQ1 MHZ(5);
#define FREQ2 MHZ(0);
#define NORUNS 5;

#define ASSERT(expr) { \
	if (!(expr)) { \
		(void) fprintf(stderr, "assertion failed (%s:%d)\n", __FILE__, __LINE__); \
		(void) abort(); \
	} \
}

/* RX is input, TX is output */
enum iodev { RX, TX };

/* common RX and TX streaming params */
struct stream_cfg {
	long long bw_hz; // Analog banwidth in Hz
	long long fs_hz; // Baseband sample rate in Hz
	long long lo_hz; // Local oscillator frequency in Hz
	const char* rfport; // Port name
};

/* static scratch mem for strings */
static char tmpstr[64];

/* IIO structs required for streaming */
static struct iio_context *ctx   = NULL;
static struct iio_channel *rx0_i = NULL;
static struct iio_channel *rx0_q = NULL;
static struct iio_channel *tx0_i = NULL;
static struct iio_channel *tx0_q = NULL;
static struct iio_buffer  *rxbuf = NULL;
static struct iio_buffer  *txbuf = NULL;

static bool stop;

/* cleanup and exit */
static void shutdown()
{
	printf("* Destroying buffers\n");
	if (rxbuf) { iio_buffer_destroy(rxbuf); }
	if (txbuf) { iio_buffer_destroy(txbuf); }

	printf("* Disabling streaming channels\n");
	if (rx0_i) { iio_channel_disable(rx0_i); }
	if (rx0_q) { iio_channel_disable(rx0_q); }
	if (tx0_i) { iio_channel_disable(tx0_i); }
	if (tx0_q) { iio_channel_disable(tx0_q); }

	printf("* Destroying context\n");
	if (ctx) { iio_context_destroy(ctx); }
	exit(0);
}

static void handle_sig(int sig)
{
	printf("Waiting for process to finish...\n");
	stop = true;
}

/* check return value of attr_write function */
static void errchk(int v, const char* what) {
	 if (v < 0) { fprintf(stderr, "Error %d writing to channel \"%s\"\nvalue may not be supported.\n", v, what); shutdown(); }
}

/* write attribute: long long int */
static void wr_ch_lli(struct iio_channel *chn, const char* what, long long val)
{
	errchk(iio_channel_attr_write_longlong(chn, what, val), what);
}

/* write attribute: string */
static void wr_ch_str(struct iio_channel *chn, const char* what, const char* str)
{
	errchk(iio_channel_attr_write(chn, what, str), what);
}

/* helper function generating channel names */
static char* get_ch_name(const char* type, int id)
{
	snprintf(tmpstr, sizeof(tmpstr), "%s%d", type, id);
	return tmpstr;
}

/* returns ad9361 phy device */
static struct iio_device* get_ad9361_phy(struct iio_context *ctx)
{
	struct iio_device *dev =  iio_context_find_device(ctx, "ad9361-phy");
	ASSERT(dev && "No ad9361-phy found");
	return dev;
}

/* finds AD9361 streaming IIO devices */
static bool get_ad9361_stream_dev(struct iio_context *ctx, enum iodev d, struct iio_device **dev)
{
	switch (d) {
	case TX: *dev = iio_context_find_device(ctx, "cf-ad9361-dds-core-lpc"); return *dev != NULL;
	case RX: *dev = iio_context_find_device(ctx, "cf-ad9361-lpc");  return *dev != NULL;
	default: ASSERT(0); return false;
	}
}

/* finds AD9361 streaming IIO channels */
static bool get_ad9361_stream_ch(struct iio_context *ctx, enum iodev d, struct iio_device *dev, int chid, struct iio_channel **chn)
{
	*chn = iio_device_find_channel(dev, get_ch_name("voltage", chid), d == TX);
	if (!*chn)
		*chn = iio_device_find_channel(dev, get_ch_name("altvoltage", chid), d == TX);
	return *chn != NULL;
}

/* finds AD9361 phy IIO configuration channel with id chid */
static bool get_phy_chan(struct iio_context *ctx, enum iodev d, int chid, struct iio_channel **chn)
{
	switch (d) {
	case RX: *chn = iio_device_find_channel(get_ad9361_phy(ctx), get_ch_name("voltage", chid), false); return *chn != NULL;
	case TX: *chn = iio_device_find_channel(get_ad9361_phy(ctx), get_ch_name("voltage", chid), true);  return *chn != NULL;
	default: ASSERT(0); return false;
	}
}

/* finds AD9361 local oscillator IIO configuration channels */
static bool get_lo_chan(struct iio_context *ctx, enum iodev d, struct iio_channel **chn)
{
	switch (d) {
	 // LO chan is always output, i.e. true
	case RX: *chn = iio_device_find_channel(get_ad9361_phy(ctx), get_ch_name("altvoltage", 0), true); return *chn != NULL;
	case TX: *chn = iio_device_find_channel(get_ad9361_phy(ctx), get_ch_name("altvoltage", 1), true); return *chn != NULL;
	default: ASSERT(0); return false;
	}
}

/* applies streaming configuration through IIO */
bool cfg_ad9361_streaming_ch(struct iio_context *ctx, struct stream_cfg *cfg, enum iodev type, int chid)
{
	struct iio_channel *chn = NULL;

	// Configure phy and lo channels
	printf("* Acquiring AD9361 phy channel %d\n", chid);
	if (!get_phy_chan(ctx, type, chid, &chn)) {	return false; }
	wr_ch_str(chn, "rf_port_select",     cfg->rfport);
	wr_ch_lli(chn, "rf_bandwidth",       cfg->bw_hz);
	wr_ch_lli(chn, "sampling_frequency", cfg->fs_hz);

	// Configure LO channel
	printf("* Acquiring AD9361 %s lo channel\n", type == TX ? "TX" : "RX");
	if (!get_lo_chan(ctx, type, &chn)) { return false; }
	wr_ch_lli(chn, "frequency", cfg->lo_hz);
	return true;
}

static ssize_t demux_sample(const struct iio_channel *chn, void *sample, size_t size, void *d){
	double val;

	iio_channel_convert(chn, &val, sample);

	return size;
}

float dither(float f)
{
	float r1 = (float)rand() / (float)RAND_MAX;
	float r2 = (float)rand() / (float)RAND_MAX;

	return f + (r1 - r2) * 16.0f;
}

void tx_thread(){
	int16_t *buf, *sine;
	int i;

	sine = malloc(sizeof(int16_t) * 1024 * 256 * 2);

	for (i = 0; i < 1024 * 256; i++) {
		sine[i*2 + 0] = dither(cos(2 * M_PI * i / 256.0) * 0x4000);
		sine[i*2 + 1] = dither(sin(2 * M_PI * i / 256.0) * 0x4000);
	}

	while (1) {
		buf = iio_buffer_start(txbuf);
		memcpy(buf, sine, 1024 * 256 * sizeof(int16_t) * 2);
		iio_buffer_push(txbuf);
	}
}

/* simple configuration and streaming */
int main (int argc, char **argv)
{
	// TX thread
	//pthread_t tx_th;
	//int thread_info;
	//void *res;
	int cnt, count;

	// File to dump data
	FILE *fp1, *fp2, *fp3;

	// Streaming devices
	struct iio_device *tx;
	struct iio_device *rx;

	// RX and TX sample counters
	size_t nrx = 0;
	size_t ntx = 0;

	// Stream configurations
	struct stream_cfg rxcfg;
	struct stream_cfg txcfg;

	ssize_t fft_size;
	fftw_complex *in, *out;
	fftw_plan plan;
	double avg, mag;
	double *out_data;

	// Listen to ctrl+c and ASSERT
	signal(SIGINT, handle_sig);

	// RX stream config
	rxcfg.bw_hz = MHZ(19.366);   	// 20 MHz rf bandwidth
	rxcfg.fs_hz = MHZ(30.72);     // 30 MS/s rx sample rate
	rxcfg.lo_hz = GHZ(1);   			// 1 GHz rf frequency
	rxcfg.rfport = "A_BALANCED"; // port A (select for rf freq.)

	// Print some device information
	printf("*RX settings\n  Bandwidth: %lld Hz\n  Baseband Sample rate: %lld Hz\n  LO frequency: %lld Hz\n", rxcfg.bw_hz, rxcfg.fs_hz, rxcfg.lo_hz);

	// TX stream config
	txcfg.bw_hz = MHZ(19.365); 	// 20 MHz rf bandwidth
	txcfg.fs_hz = MHZ(30.72);   // 30 MS/s tx sample rate
	txcfg.lo_hz = GHZ(1); 			// 1 GHz rf frequency
	txcfg.rfport = "A"; // port A (select for rf freq.)

	printf("* Acquiring IIO context\n");
	//ASSERT((ctx = iio_create_default_context()) && "No context");
	ASSERT((ctx = iio_create_network_context("192.168.1.227")) && "No context");
	ASSERT(iio_context_get_devices_count(ctx) > 0 && "No devices");

	printf("* Acquiring AD9361 streaming devices\n");
	ASSERT(get_ad9361_stream_dev(ctx, TX, &tx) && "No tx dev found");
	ASSERT(get_ad9361_stream_dev(ctx, RX, &rx) && "No rx dev found");

	printf("* Configuring AD9361 for streaming\n");
	ASSERT(cfg_ad9361_streaming_ch(ctx, &rxcfg, RX, 0) && "RX port 0 not found");
	ASSERT(cfg_ad9361_streaming_ch(ctx, &txcfg, TX, 0) && "TX port 0 not found");

	printf("* Initializing AD9361 IIO streaming channels\n");
	ASSERT(get_ad9361_stream_ch(ctx, RX, rx, 0, &rx0_i) && "RX chan i not found");
	ASSERT(get_ad9361_stream_ch(ctx, RX, rx, 1, &rx0_q) && "RX chan q not found");
	ASSERT(get_ad9361_stream_ch(ctx, TX, tx, 0, &tx0_i) && "TX chan i not found");
	ASSERT(get_ad9361_stream_ch(ctx, TX, tx, 1, &tx0_q) && "TX chan q not found");

	printf("* Number of RX channels: %d\n", iio_device_get_channels_count(rx));

	printf("* Enabling IIO streaming channels\n");
	iio_channel_enable(rx0_i);
	iio_channel_enable(rx0_q);
	iio_channel_enable(tx0_i);
	iio_channel_enable(tx0_q);

	int buffer_size = 1024 * 1024;

	printf("* Creating non-cyclic IIO buffers with 1 MiS\n");
	rxbuf = iio_device_create_buffer(rx, buffer_size, false);
	if (!rxbuf) {
		perror("Could not create RX buffer");
		shutdown();
	}
	txbuf = iio_device_create_buffer(tx, buffer_size, false);
	if (!txbuf) {
		perror("Could not create TX buffer");
		shutdown();
	}

	// configure fft
  fft_size = 16384;	// 32768; Same size as iio_buffer
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fft_size);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fft_size);
	out_data = malloc(sizeof(double)*fft_size);
	plan = fftw_plan_dft_1d(fft_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

	printf("* Starting IO streaming (press CTRL+C to cancel)\n");

	fp1 = fopen("output.csv", "w+");
	fp2 = fopen("input.csv", "w+");
	// Create TX thread
	//pthread_create (&tx_th, NULL, (void*) &tx_thread, NULL);
	count = NORUNS;

	while (!stop && count > 0){
		ssize_t nbytes_rx, nbytes_tx;
		char *p_dat, *p_end;
		ptrdiff_t p_inc;

		// Schedule TX buffer
		nbytes_tx = iio_buffer_push(txbuf);
		if (nbytes_tx < 0) { printf("Error pushing buf %d\n", (int) nbytes_tx); shutdown(); }

		// Refill RX buffer
		nbytes_rx = iio_buffer_refill(rxbuf);
		if (nbytes_rx < 0) { printf("Error refilling buf %d\n",(int) nbytes_rx); shutdown(); }

		// READ: Get pointers to RX buf and read IQ from RX buf port 0

		p_inc = iio_buffer_step(rxbuf);
		p_end = iio_buffer_end(rxbuf);

		// Convert to native format
		iio_buffer_foreach_sample(rxbuf, demux_sample, NULL);

		// Dump received data to file for analysis
		cnt = 0;
		for (p_dat = (char *)iio_buffer_first(rxbuf, rx0_i); p_dat < p_end; p_dat += p_inc) {
			// Get I and Q and save to file
			const int16_t i = ((int16_t*)p_dat)[0]; // Real (I)
			const int16_t q = ((int16_t*)p_dat)[1]; // Imag (Q)

			// Copy captured data into fftw3 in buffer
			if ( cnt < fft_size ){
				in[cnt] = i + q*I;
				cnt++;
			}
			// Print data to file
			fprintf(fp2, "%d,%d\n", i, q);
		}

		fftw_execute(plan);
		avg = AVERAGES;
		if (avg && avg != 128 )
			avg = 1.0f / avg;

		// Sample counter increment and status output
		nrx += nbytes_rx / iio_device_get_sample_size(rx);
		ntx += nbytes_tx / iio_device_get_sample_size(tx);
		printf("\tRX %8.2f MSmp, TX %8.2f MSmp\n", nrx/1e6, ntx/1e6);

		fp3 = fopen("fft.csv", "w");
		for(cnt = 0; cnt<fft_size; cnt++){
			mag = 10*log10( (creal(out[cnt]) * creal(out[cnt]) + cimag(out[cnt]) * cimag(out[cnt])) / ((unsigned long long)fft_size * fft_size));

			if (out_data[cnt] == FLT_MAX) {
				/* Don't average the first iteration */
				 out_data[cnt] = mag;
			} else if (!avg) {
				/* keep peaks */
				if (out_data[cnt] <= mag)
					out_data[cnt] = mag;
			} else if (avg == 128) {
				/* keep min */
				if (out_data[cnt] >= mag)
					out_data[cnt] = mag;
			} else {
				/* do an average */
				out_data[cnt] = ((1 - avg) * out_data[cnt]) + (avg * mag);
			}
			//fprintf(fp3, "%lf,%lf\n", out[0][cnt], out[1][cnt]);
			//fprintf(fp3, "%lf\n", 10*log10( (out[0][cnt]*out[0][cnt] + out[1][cnt]*out[1][cnt]) / (fft_size*fft_size) ) );
			fprintf(fp3, "%lf\n", out_data[cnt]);

		}
		fclose(fp3);

		// WRITE: Get pointers to TX buf and write IQ to TX buf port 0
		p_inc = iio_buffer_step(txbuf);
		p_end = iio_buffer_end(txbuf);

		float freq1 = 2.0 * M_PI * FREQ1; // sine wave (should be easy to spot);
		float freq2 = 2.0 * M_PI * FREQ2; // sine wave (should be easy to spot);

		double ampl = 32767;

		double i = 1. / txcfg.fs_hz;

		for (p_dat = iio_buffer_first(txbuf, tx0_i); p_dat < p_end; p_dat += p_inc) {
			// fill with sine wave
			short ipart = ampl * sin(freq1 * i);// + ampl * sin(freq2 * i);
			short qpart = ampl * cos(freq1 * i);// + ampl * cos(freq2 * i);

			((int16_t *)p_dat)[0] = ipart;
			((int16_t *)p_dat)[1] = qpart;

			// Save what's actually in the TX buffer to file
			fprintf(fp1, "%d,%d\n", ((int16_t*)p_dat)[0], ((int16_t*)p_dat)[1]);

			i += 1. / txcfg.fs_hz;
		}
		count--;
	}

	// thread_info = pthread_cancel(tx_th);
  // if (thread_info != 0)
  // 	printf("pthread_cancel error\n");
	//
  // /* Join with thread to see what its exit status was */
	// thread_info = pthread_join(tx_th, &res);
  // if (thread_info != 0)
  // 	printf("pthread_join error\n");
	fclose(fp1);
	fclose(fp2);
	fftw_destroy_plan(plan);
	fftw_free(in);
	fftw_free(out);

	// Temp, quit now as hing on buffer destroy? Need to figure out why. mem leakage :-/
	//return (0);
	shutdown();

	return 0;
}
