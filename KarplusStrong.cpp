// Karl Yerkes
// 2023-01-30
// MAT 240B ~ Audio Programming
// Assignment 3 ~ Karplus-Strong string modeling
//

#include <juce_audio_processors/juce_audio_processors.h>

template <typename T>
T mtof(T m) {
  return T(440) * pow(T(2), (m - T(69)) / T(12));
}
template <typename T>
T dbtoa(T db) {
  return pow(T(10), db / T(20));
}

// valid on (-1, 1)
template <class T>
inline T sine(T n) {
  T nn = n * n;
  return n * (T(3.138982) +
              nn * (T(-5.133625) + nn * (T(2.428288) - nn * T(0.433645))));
}

template <class T>
inline T softclip(T x) {
  if (x >= T(1)) return T(1);
  if (x <= T(-1)) return T(-1);
  return (T(3) * x - x * x * x) / T(2);
}

template <class T>
inline T wrap(T v, T hi, T lo) {
  if (lo == hi) return lo;

  // if(v >= hi){
  if (!(v < hi)) {
    T diff = hi - lo;
    v -= diff;
    if (!(v < hi)) v -= diff * (T)(unsigned)((v - lo) / diff);
  } else if (v < lo) {
    T diff = hi - lo;
    v += diff;  // this might give diff if range is too large, so check at end
                // of block...
    if (v < lo) v += diff * (T)(unsigned)(((lo - v) / diff) + 1);
    if (v == diff) return std::nextafter(v, lo);
  }
  return v;
}

struct BooleanOscillator {
  float value = 1;
  float increment = 0;
  void frequency(float hertz, float samplerate) {
    assert(hertz >= 0);
    increment = hertz / samplerate;
  }
  void period(float hertz, float samplerate) {
    frequency(1 / hertz, samplerate);
  }
  bool operator()() {
    value += increment;
    bool b = value >= 1;
    value = wrap(value, 1.0f, 0.0f);
    return b;
  }
};

// https://en.wikipedia.org/wiki/Harmonic_oscillator
struct MassSpringModel {
  // this the whole state of the simulation
  //
  float position{0};  // m
  float velocity{0};  // m/s

  // These are cached properties of the model; They govern the behaviour. We
  // recalculate them given frequency, decay time, and playback rate.
  //
  float springConstant{0};      // N/m
  float dampingCoefficient{0};  // N·s/m

  void show() {
    printf("position:%f velocity:%f springConstant:%f dampingCoefficient:%f\n",
           position, velocity, springConstant, dampingCoefficient);
  }
  void reset() {
    // show();
    position = velocity = 0;
  }

  float operator()() {
    // This is semi-implicit Euler integration with time-step 1. The
    // playback rate is "baked into" the constants. Spring force and damping
    // force are accumulated into velocity. We let mass is 1, so it
    // disappears. Velocity is accumulated into position which is
    // interpreted as oscillator amplitude.
    //
    float acceleration = 0;

    // XXX put code here

    velocity += acceleration;
    position += velocity;

    /*
        printf("position:%f velocity:%f springConstant:%f
       dampingCoefficient: %f\n", position, velocity, springConstant,
       dampingCoefficient);
    */
    return position;
  }

  // Use these to measure the kinetic, potential, and total energy of the
  // system.
  float ke() { return velocity * velocity / 2; }
  float pe() { return position * position * springConstant / 2; }
  float te() { return ke() + pe(); }

  // "Kick" the mass-spring system such that we get a nice (-1, 1) oscillation.
  //
  void trigger() {
    // We want the "mass" to move in (-1, 1). What is the potential energy
    // of a mass-spring system at 1? PE == k * x * x / 2 == k / 2. So, we
    // want a system with k / 2 energy, but we don't want to just set the
    // displacement to 1 because that would make a click. Instead, we want
    // to set the velocity. What velocity would we need to have energy k /
    // 2? KE == m * v * v / 2 == k / 2. or v * v == k. so...
    //

    // XXX put code here

    // How might we improve on this? Consider triggering at a level
    // depending on frequency according to the Fletcher-Munson curves.
  }

  void recalculate(float frequency, float decayTime, float playbackRate) {
    // sample rate is "baked into" these constants to save on per-sample
    // operations.

    // XXX put code here
  }
};

class DelayLine : public std::vector<float> {
  //
  int index = 0;

 public:
  float read(float samples_ago) {
    jassert(samples_ago < size());
    float i = index - samples_ago;
    if (i < 0) {
      i += size();
    }
    return at((int)i);  // no linear interpolation
  }

  void write(float value) {
    jassert(size() > 0);
    at(index) = value;  // overwrite the oldest value

    // handle the wrapping for circular buffer
    index++;
    if (index >= size()) index = 0;
  }

  void allocate(float seconds, float samplerate) {
    // floor(seconds * samplerate) + 1 samples
    resize((int)floor(seconds * samplerate) + 1);
  }
};

class BiquadFilter {
  // Audio EQ Cookbook
  // http://www.musicdsp.org/files/Audio-EQ-Cookbook.txt

  // x[n-1], x[n-2], y[n-1], y[n-2]
  float x1 = 0, x2 = 0, y1 = 0, y2 = 0;

  // filter coefficients
  float b0 = 1, b1 = 0, b2 = 0, a1 = 0, a2 = 0;

 public:
  float operator()(float x0) {
    // Direct Form 1, normalized...
    float y0 = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
    y2 = y1;
    y1 = y0;
    x2 = x1;
    x1 = x0;
    return y0;
  }

  void normalize(float a0) {
    b0 /= a0;
    b1 /= a0;
    b2 /= a0;
    a1 /= a0;
    a2 /= a0;
  }

  void lpf(float f0, float Q, float samplerate) {
    float w0 = 2 * float(M_PI) * f0 / samplerate;
    float alpha = sin(w0) / (2 * Q);
    b0 = (1 - cos(w0)) / 2;
    b1 = 1 - cos(w0);
    b2 = (1 - cos(w0)) / 2;
    float a0 = 1 + alpha;
    a1 = -2 * cos(w0);
    a2 = 1 - alpha;

    normalize(a0);
  }
};

struct KarplusStrongModel {
  DelayLine delay;
  BiquadFilter filter;
  float delayTime = 0;  // in *samples*
  float feedback = 1;

  void init(float seconds, float samplerate) {
    delay.allocate(seconds, samplerate);
    zero();
  }

  void configure(float hertz, float seconds, float samplerate) {
    delayTime = samplerate / hertz;
    filter.lpf(1000, 0.4, samplerate);
  }

  float trigger() {
    static int random = 0;
    for (int i = 0; i < delay.size(); ++i) {
      random += 12345;
      random *= 1103515245;
      delay[i] = random / 2147483647.0f;
    }
  }

  float zero() {
    for (int i = 0; i < delay.size(); ++i) {
      delay[i] = 0;
    }
  }

  float operator()() {
    float v = filter(delay.read(delayTime));
    delay.write(v * feedback);
    return v;
  }
};

using namespace juce;

class KarplusStrong : public AudioProcessor {
  AudioParameterFloat* gain;
  AudioParameterFloat* note;
  AudioParameterFloat* rate;
  AudioParameterFloat* cutoff;
  AudioParameterFloat* resonance;
  AudioParameterFloat* feedback;

  BooleanOscillator timer;
  // MassSpringModel string;
  KarplusStrongModel string;

  /// add parameters here ///////////////////////////////////////////////////

 public:
  KarplusStrong()
      : AudioProcessor(BusesProperties()
                           .withInput("Input", AudioChannelSet::stereo())
                           .withOutput("Output", AudioChannelSet::stereo())) {
    addParameter(gain = new AudioParameterFloat(
                     {"gain", 1}, "Gain",
                     NormalisableRange<float>(-65, -1, 0.01f), -65));
    addParameter(
        note = new AudioParameterFloat(
            {"note", 1}, "Note", NormalisableRange<float>(-2, 129, 0.01f), 40));
    addParameter(
        rate = new AudioParameterFloat(
            {"rate", 1}, "Rate", NormalisableRange<float>(-40, 40, 0.01f), -4));
    addParameter(cutoff = new AudioParameterFloat(
                     {"cutoff", 1}, "Cutoff",
                     NormalisableRange<float>(-2, 129, 0.01f), 100));
    addParameter(resonance = new AudioParameterFloat(
                     {"resonance", 1}, "Resonance",
                     NormalisableRange<float>(0.00001, 4, 0.01f), 0.6));
    addParameter(feedback = new AudioParameterFloat(
                     {"feedback", 1}, "Feedback",
                     NormalisableRange<float>(0.0, 2.0, 0.01f), 1.0));

    /// add parameters here /////////////////////////////////////////////

    // XXX juce::getSampleRate() is not valid here
  }

  float previous = 0;

  /// handling the actual audio! ////////////////////////////////////////////
  void processBlock(AudioBuffer<float>& buffer, MidiBuffer&) override {
    buffer.clear(0, 0, buffer.getNumSamples());
    auto left = buffer.getWritePointer(0, 0);
    auto right = buffer.getWritePointer(1, 0);

    timer.frequency(mtof(rate->get()), getSampleRate());
    string.feedback = feedback->get();
    string.filter.lpf(mtof(cutoff->get()), resonance->get(), getSampleRate());
    for (int i = 0; i < buffer.getNumSamples(); ++i) {
      if (timer()) {
        string.configure(mtof(note->get()), 1, getSampleRate());
        string.trigger();
      }

      left[i] = previous = string() * dbtoa(gain->get());
      right[i] = left[i];
    }
  }

  /// handle doubles ? //////////////////////////////////////////////////////
  // void processBlock(AudioBuffer<double>& buffer, MidiBuffer&) override {
  //   buffer.applyGain(dbtoa((float)*gain));
  // }

  /// start and shutdown callbacks///////////////////////////////////////////
  void prepareToPlay(double samplerate, int) override {
    // XXX when does this get called? seems to not get called in stand-alone
    string.init(1, samplerate);
  }
  void releaseResources() override {}

  /// maintaining persistant state on suspend ///////////////////////////////
  void getStateInformation(MemoryBlock& destData) override {
    MemoryOutputStream(destData, true).writeFloat(*gain);
  }

  void setStateInformation(const void* data, int sizeInBytes) override {
    gain->setValueNotifyingHost(
        MemoryInputStream(data, static_cast<size_t>(sizeInBytes), false)
            .readFloat());
  }

  /// general configuration /////////////////////////////////////////////////
  const String getName() const override { return "Quasi Band Limited"; }
  double getTailLengthSeconds() const override { return 0; }
  bool acceptsMidi() const override { return true; }
  bool producesMidi() const override { return false; }

  /// for handling presets //////////////////////////////////////////////////
  int getNumPrograms() override { return 1; }
  int getCurrentProgram() override { return 0; }
  void setCurrentProgram(int) override {}
  const String getProgramName(int) override { return "None"; }
  void changeProgramName(int, const String&) override {}

  /// ?????? ////////////////////////////////////////////////////////////////
  bool isBusesLayoutSupported(const BusesLayout& layouts) const override {
    const auto& mainInLayout = layouts.getChannelSet(true, 0);
    const auto& mainOutLayout = layouts.getChannelSet(false, 0);

    return (mainInLayout == mainOutLayout && (!mainInLayout.isDisabled()));
  }

  /// automagic user interface //////////////////////////////////////////////
  AudioProcessorEditor* createEditor() override {
    return new GenericAudioProcessorEditor(*this);
  }
  bool hasEditor() const override { return true; }

 private:
  JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(KarplusStrong)
};

AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
  return new KarplusStrong();
}