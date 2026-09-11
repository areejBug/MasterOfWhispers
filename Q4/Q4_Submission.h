#ifndef Q4_SUBMISSION_H
#define Q4_SUBMISSION_H

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

#pragma pack(push,1)

struct WavHeader {
char chunkId[4];
int chunkSize;
char format[4];
char subchunk1Id[4];
int subchunk1Size;
short audioFormat;
short numChannels;
int sampleRate;
int byteRate;
short blockAlign;
short bitsPerSample;
char subchunk2Id[4];
int subchunk2Size;
};

#pragma pack(pop)

class Complex {
public:
    double re, im;
    Complex(double r = 0, double i = 0);
    Complex operator+(const Complex&)const;
    Complex operator-(const Complex&)const;
    Complex operator*(const Complex&)const;
};

class Filter {
public:
    float* w;int n;
    Filter(int n);
    Filter(const Filter& other);
    Filter& operator=(const Filter& other);
    ~Filter();
    static Filter low(int n, int sr, double c);
    static Filter high(int n, int sr, double c);
    static Filter band(int n, int sr, double l, double h);
};

class Spectrum;

class Whisper {
public:
    WavHeader header;
    short* samples;
    int n;

    Whisper();
    Whisper(const char*);
    ~Whisper();

    void normalize();
    Whisper& operator*=(double);
    void reverse();
    Whisper operator+(const Whisper&);
    Whisper& operator+=(const Whisper&);
    Whisper splice(double startSec, double endSec);
    Spectrum FFT();
    void save(const char*);

    Whisper(const Whisper& other) {
        header = other.header;
        n = other.n;
        samples = new short[n];
        for (int i = 0;i < n;i++) {
            samples[i] = other.samples[i];
        }
    }

    Whisper& operator=(const Whisper& other) {
        if (this != &other) {
            delete[] samples;
            header = other.header;
            n = other.n;
            samples = new short[n];
            for (int i = 0;i < n;i++) {
                samples[i] = other.samples[i];
            }
        }
        return *this;
    }
};

class Spectrum {
public:
    Complex* data;int n, orig, rate;
    WavHeader header;
    Spectrum(int);
    Spectrum(const Spectrum& other);
    Spectrum& operator=(const Spectrum& other);
    ~Spectrum();
    void apply(Filter&);
    Whisper IFFT();
};

#endif
