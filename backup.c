/*
 *
 * AUTHOR: David Scott (On behalf of IQHQ)
 * PROGRAM: Spectrum analyser for AD9361 using libiio
 *
 * NOTES:
 * The frequency resolution f_res,  is equal to the reciprocal of your time window duration T_win.
 * So if you have M samples going into your FFT and your sampling period is T seconds, then T_win = M*T,
 * thus f_res = 1/(M*T). So the frequency of the kth bin (in Hz), where k = 0 corresponds to dc,
 * will be k/(M*T) or (k/M)*F_s, if F_s is your sampling frequency, with F_s = 1/T.
 *
 * sampling period T = 1 / 30.72e6 = 0.000000033 = 33e-9 (s) i.e. 33ns
 * number of samples M = 1024*1024 = 1048576 samples
 * > T_win = M*T = 1048576*33e-9 = 0.034603008 (s)
 * > f_res = 1/0.034603008 = 28.899221709
 *
 * i.e freq at any given bin = (k/M)*F_s
 *
*/

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
#include <time.h>

#ifdef __APPLE__
#include <iio/iio.h>
#else
#include <iio.h>
#endif

/* FFT settings
	FFT timing based options
	* FFTW_ESTIMATE 	- Fastest but least accurate
	* FFTW_MEASURE 		- slower but more accurate
	* FFTW_PATIENT 		- even slower but even more accurate
	* FFTW_EXHAUSTIVE - slowest but most accurate
*/
#define AVERAGES 100
#define FFT_SIZE 1024
#define FFTW_TYPE FFTW_MEASURE

/* Buffer settings */
#define IN_BUFFER_SIZE FFT_SIZE * AVERAGES // Read enough to fill all the fft average buffers
#define OUT_BUFFER_SIZE 1024

/* Test signal settings */
#define FREQ1 MHZ(1)
#define FREQ2 MHZ(4)
#define FREQ3 MHZ(7)
#define AMP1 384
#define AMP2 384
#define AMP3 384

/* RX Spectrum settings */
#define START_SWEEP GHZ(1)
#define END_SWEEP GHZ(2)
#define RX_BW MHZ(20)
#define RX_FS MHZ(30.72)
#define RX_LO GHZ(2) // This will get changed when sweeping

/* RX Spectrum settings */
#define TX_BW MHZ(20)
#define TX_FS MHZ(30.72)
#define TX_LO GHZ(2) // This will get changed for sweeping

/* General settings */
#define NETWORK_DEVICE "192.168.1.227"
#define RUN_TIME 2

/* Helper macros */
#define KHZ(x) ((long long)(x*1000.0 + .5))
#define MHZ(x) ((long long)(x*1000000.0 + .5))
#define GHZ(x) ((long long)(x*1000000000.0 + .5))

/* Assertation helper */
#define ASSERT(expr) { \
	if (!(expr)) { \
		(void) fprintf(stderr, "assertion failed (%s:%d)\n", __FILE__, __LINE__); \
		(void) abort(); \
	} \
}

/* RX is input, TX is output */
enum iodev { RX, TX };

/* Common RX and TX streaming params */
struct stream_cfg {
	long long bw_hz; 		// Analog banwidth in Hz
	long long fs_hz; 		// Baseband sample rate in Hz
	long long lo_hz; 		// Local oscillator frequency in Hz
	const char* rfport; // Port name
};

/* Static scratch mem for strings */
static char tmpstr[64];

/* IIO structs required for streaming */
static struct iio_context *ctx   = NULL;
static struct iio_channel *rx0_i = NULL;
static struct iio_channel *rx0_q = NULL;
static struct iio_channel *tx0_i = NULL;
static struct iio_channel *tx0_q = NULL;
static struct iio_buffer  *rxbuf = NULL;
static struct iio_buffer  *txbuf = NULL;

/* global variables */
/* Streaming devices */
struct iio_device *tx;
struct iio_device *rx;
/* RX and TX sample counters */
size_t nrx = 0;
size_t ntx = 0;
/* Stream configurations */
struct stream_cfg rxcfg;
struct stream_cfg txcfg;

/* FFT variables */
ssize_t fft_size;
fftw_complex *avg_buffers;
fftw_complex *in, *out;
fftw_plan plan;
double *win;
fftw_complex total, average;

/* General global variables */
FILE *fp1, *fp2, *fp3;
static bool stop;
int nbytes_rx;
clock_t start, end;
double cpu_time_used, avg_cpu_time_used;

/* Cleanup and exit */
static void shutdown()
{
	printf("* Freeing file pointers\n");
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);

	printf("* Destroying buffers\n");
	fftw_destroy_plan(plan);
	fftw_free(in);
	fftw_free(out);
	fftw_free(avg_buffers);
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

/* Check return value of attr_write function */
static void errchk(int v, const char* what) {
	 if (v < 0) { fprintf(stderr, "Error %d writing to channel \"%s\"\nvalue may not be supported.\n", v, what); shutdown(); }
}

/* Write attribute: long long int */
static void wr_ch_lli(struct iio_channel *chn, const char* what, long long val)
{
	errchk(iio_channel_attr_write_longlong(chn, what, val), what);
}

/* Write attribute: string */
static void wr_ch_str(struct iio_channel *chn, const char* what, const char* str)
{
	errchk(iio_channel_attr_write(chn, what, str), what);
}

/* Helper function generating channel names */
static char* get_ch_name(const char* type, int id)
{
	snprintf(tmpstr, sizeof(tmpstr), "%s%d", type, id);
	return tmpstr;
}

/* Returns ad9361 phy device */
static struct iio_device* get_ad9361_phy(struct iio_context *ctx)
{
	struct iio_device *dev =  iio_context_find_device(ctx, "ad9361-phy");
	ASSERT(dev && "No ad9361-phy found");
	return dev;
}

/* Finds AD9361 streaming IIO devices */
static bool get_ad9361_stream_dev(struct iio_context *ctx, enum iodev d, struct iio_device **dev)
{
	switch (d) {
	case TX: *dev = iio_context_find_device(ctx, "cf-ad9361-dds-core-lpc"); return *dev != NULL;
	case RX: *dev = iio_context_find_device(ctx, "cf-ad9361-lpc");  return *dev != NULL;
	default: ASSERT(0); return false;
	}
}

/* Finds AD9361 streaming IIO channels */
static bool get_ad9361_stream_ch(struct iio_context *ctx, enum iodev d, struct iio_device *dev, int chid, struct iio_channel **chn)
{
	*chn = iio_device_find_channel(dev, get_ch_name("voltage", chid), d == TX);
	if (!*chn)
		*chn = iio_device_find_channel(dev, get_ch_name("altvoltage", chid), d == TX);
	return *chn != NULL;
}

/* Finds AD9361 phy IIO configuration channel with id chid */
static bool get_phy_chan(struct iio_context *ctx, enum iodev d, int chid, struct iio_channel **chn)
{
	switch (d) {
	case RX: *chn = iio_device_find_channel(get_ad9361_phy(ctx), get_ch_name("voltage", chid), false); return *chn != NULL;
	case TX: *chn = iio_device_find_channel(get_ad9361_phy(ctx), get_ch_name("voltage", chid), true);  return *chn != NULL;
	default: ASSERT(0); return false;
	}
}

/* Finds AD9361 local oscillator IIO configuration channels */
static bool get_lo_chan(struct iio_context *ctx, enum iodev d, struct iio_channel **chn)
{
	switch (d) {
	 // LO chan is always output, i.e. true
	case RX: *chn = iio_device_find_channel(get_ad9361_phy(ctx), get_ch_name("altvoltage", 0), true); return *chn != NULL;
	case TX: *chn = iio_device_find_channel(get_ad9361_phy(ctx), get_ch_name("altvoltage", 1), true); return *chn != NULL;
	default: ASSERT(0); return false;
	}
}

/* Applies streaming configuration through IIO */
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

// Generte Hanning window
static double win_hanning(int j, int n)
{
	double a = 2.0 * M_PI / (n - 1), w;
	w = 0.5 * (1.0 - cos(a * j));
	return (w);
}

// Seperate thread for receive chain
void rx_thread(){
	// General variables
	int cnt, avg_cnt;

	while (!stop){
		char *p_dat, *p_end;
		ptrdiff_t p_inc;

		// Refill RX buffer
		nbytes_rx = iio_buffer_refill(rxbuf);
		if (nbytes_rx < 0) { printf("Error refilling buf %d\n",(int) nbytes_rx); shutdown(); }

		// READ: Get pointers to RX buf and read IQ from RX buf port 0
		p_inc = iio_buffer_step(rxbuf);
		p_end = iio_buffer_end(rxbuf);

		cnt = 0;
		for (p_dat = (char *)iio_buffer_first(rxbuf, rx0_i); p_dat < p_end; p_dat += p_inc) {
			// Get I and Q and save to file
			const int16_t i = ((int16_t*)p_dat)[0]; // Real (I)
			const int16_t q = ((int16_t*)p_dat)[1]; // Imag (Q)

			//uncomment next line to save received data to input.csv file
			//fprintf(fp2, "%d,%d\n", ((int16_t*)p_dat)[0], ((int16_t*)p_dat)[1]);

			// copies AVERAGES * fft_size data into buffer
			if ( cnt < fft_size*AVERAGES ){
				avg_buffers[cnt] = i + (I*q);
				cnt++;
			}
		}

		// place fft_size samples into fftw3 in buffer, call FFT and save results for averaging
		for ( avg_cnt = 0; avg_cnt < AVERAGES; avg_cnt++ ){
			for ( cnt = avg_cnt*fft_size; cnt < (avg_cnt*fft_size)+fft_size; cnt++ ){
				in[cnt-(avg_cnt*fft_size)] = avg_buffers[cnt]*win[cnt-(avg_cnt*fft_size)];
			}

			start = clock();
			fftw_execute(plan);
			end = clock();
			cpu_time_used += ((double)(end-start)) / CLOCKS_PER_SEC;

			// copy data back into avg buffer
			for ( cnt = 0; cnt < fft_size; cnt++ ){
				avg_buffers[(avg_cnt*fft_size)+cnt] = out[cnt];
			}
		}

		avg_cpu_time_used = cpu_time_used/AVERAGES;

		// Average amplitude of results
		for ( cnt = 0; cnt < fft_size; cnt++ ){
			total = 0;
			for ( avg_cnt = 0; avg_cnt < AVERAGES; avg_cnt++ ){
				total += cabs(avg_buffers[(avg_cnt*fft_size)+cnt]);
			}
			average = total/AVERAGES;
			out[cnt] = average;
		}

		// Discard first couple of bins, just to limit 1/f noise (not sure if this should be done or not!)
		out[0] = 0;
		out[1] = 0;

		rewind(fp3);
		for(cnt = 0; cnt<fft_size; cnt++){
			fprintf(fp3, "%lf %lf\n", (((double)cnt/(double)fft_size)*(double)rxcfg.fs_hz)/1e6, 20*log10((cabs(out[cnt]))/(fft_size*fft_size) ));
			//fprintf(fp, "%lf %lf\n", (((double)cnt/(double)fft_size)*(double)rxcfg.fs_hz), cabs(out[cnt]) );
			//fprintf(fp, "%lf %lf\n", (((double)cnt/(double)fft_size)*(double)rxcfg.fs_hz), 10*log10(creal(out[cnt])*creal(out[cnt]) + cimag(out[cnt])*cimag(out[cnt])) / ((unsigned long long)fft_size*fft_size) );
		}
	}
}

/*
Set up transmission of test sine tones. Cyclic buffers mean that
this is a one shot deal. Works MUCH better than constantly pushing
an output buffer which essentially holds the same data each time!
currently supports up to 3 sine tones as test signals
*/
void test_signals(float fr1, float fr2, float fr3, double ampl1, double ampl2, double ampl3){
	char *p_dat, *p_end;
	ptrdiff_t p_inc;
	int nbytes_tx;

	float freq1 = 2.0 * M_PI * fr1;
	float freq2 = 2.0 * M_PI * fr2;
	float freq3 = 2.0 * M_PI * fr3;

	double i;
	short qpart, ipart;

	// WRITE: Get pointers to TX buf and write IQ to TX buf port 0
	p_inc = iio_buffer_step(txbuf);
	p_end = iio_buffer_end(txbuf);

	i = 1. / txcfg.fs_hz;

	for (p_dat = iio_buffer_first(txbuf, tx0_i); p_dat < p_end; p_dat += p_inc) {
	// fill with sine wave
	ipart = (ampl1 * cos(freq1 * i)) + (ampl2 * cos(freq2 * i)) + (ampl3 * cos(freq3 * i));
	qpart = (ampl1 * sin(freq1 * i)) + (ampl2 * sin(freq2 * i)) + (ampl3 * sin(freq3 * i));

	((int16_t *)p_dat)[0] = ipart;
	((int16_t *)p_dat)[1] = qpart;
	i += 1. / txcfg.fs_hz;

	//Uncomment next line to save what's actually in the TX buffer to file
	//fprintf(fp1, "%d,%d\n", ((int16_t*)p_dat)[0], ((int16_t*)p_dat)[1]);
	}

	// Schedule TX buffer
	nbytes_tx = iio_buffer_push(txbuf);
	if (nbytes_tx < 0) { printf("Error pushing buf %d\n", (int) nbytes_tx); shutdown(); }
}

void setup_device(){
	int buffer_size;

	// RX stream config
	rxcfg.bw_hz = RX_BW;
	rxcfg.fs_hz = RX_FS;
	rxcfg.lo_hz = RX_LO;
	rxcfg.rfport = "A_BALANCED"; // port A (select for rf freq.)

	// TX stream config
	txcfg.bw_hz = TX_BW;
	txcfg.fs_hz = TX_FS;
	txcfg.lo_hz = TX_LO;
	txcfg.rfport = "A"; // port A (select for rf freq.)

	printf("* Acquiring IIO context\n");
	ASSERT((ctx = iio_create_network_context(NETWORK_DEVICE)) && "No context");
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

	buffer_size = IN_BUFFER_SIZE;

	printf("* Creating non-cyclic IIO buffers with 1 MiS\n");
	rxbuf = iio_device_create_buffer(rx, buffer_size, false);
	if (!rxbuf) {
		perror("Could not create RX buffer");
		shutdown();
	}

	buffer_size = OUT_BUFFER_SIZE;
	txbuf = iio_device_create_buffer(tx, buffer_size, true);
	if (!txbuf) {
		perror("Could not create TX buffer");
		shutdown();
	}
}

/* Main entry point */
int main (int argc, char **argv)
{
	// General variables
	int cnt;

	// TX thread
	pthread_t rx_th;
	int thread_info, timer;
	void *res;

	// Listen to ctrl+c and ASSERT
	signal(SIGINT, handle_sig);

	// Configure AD9361 Device
	setup_device();

	// configure fft
  fft_size = FFT_SIZE; //16384;	// Same size as iio_buffer
	in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fft_size);
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fft_size);
	avg_buffers = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*fft_size*AVERAGES);
	plan = fftw_plan_dft_1d(fft_size, in, out, FFTW_FORWARD, FFTW_TYPE);
	win = malloc(sizeof(double)*fft_size);

	// Define Hanning window
	for (cnt = 0; cnt < fft_size; cnt ++)
		win[cnt] = win_hanning(cnt, fft_size);

	printf("* Starting IO streaming (press CTRL+C to cancel)\n");
	fp1 = fopen("output.csv", "w+");
	fp2 = fopen("input.csv", "w+");
	fp3 = fopen("fft.csv", "w");

	// Call TX function
	test_signals(FREQ1, FREQ2, FREQ3, AMP1, AMP2, AMP3);

	// Create RX threads
	pthread_create (&rx_th, NULL, (void*) &rx_thread, NULL);

	// Loop for a times
	for ( timer = 0; timer < RUN_TIME; timer++ ){
		sleep(1);
		nrx += nbytes_rx / iio_device_get_sample_size(rx);
		printf("\tRX %8.2f MSmp\n", nrx/1e6);
	}

	// call cancel flag for threads
	stop = true;

  // Join with thread to see what its exit status was
	thread_info = pthread_join(rx_th, &res);
	if (thread_info != 0)
	 	printf("pthread_join error\n");

	printf("Average time required for FFT  = %.2lf (us)\n", avg_cpu_time_used*1e6);
	printf("* starting shutdown...\n" );
	shutdown();

	return 0;
}
