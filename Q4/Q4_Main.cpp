#include "Q4_Submission.h"

const double PI = 3.141592653589793;
const short int S16_MAX = 32767;
const short int S16_MIN = -32768;

int nextPow2(int n) { 
    int p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}
short int clip(double v) {
    if (v > S16_MAX)
        return S16_MAX;
    if (v < S16_MIN)
        return S16_MIN;
    return (short int)v;
}

Complex::Complex(double r, double i) :re(r), im(i) {}
Complex Complex::operator+(const Complex& o)const {
    return Complex(re + o.re, im + o.im); 
}
Complex Complex::operator-(const Complex& o)const { 
    return Complex(re - o.re, im - o.im); 
}
Complex Complex::operator*(const Complex& o)const { 
    return Complex(re * o.re - im * o.im, re * o.im + im * o.re);
}

Filter::Filter(int n) :n(n) {
    w = new float[n];
}
Filter::Filter(const Filter& other) : n(other.n) {
    w = new float[n];
    for (int i = 0;i < n;i++) {
        w[i] = other.w[i];
    }
}
Filter& Filter::operator=(const Filter& other) {
    if (this != &other) {
        delete[] w;
        n = other.n;
        w = new float[n];
        for (int i = 0;i < n;i++) {
            w[i] = other.w[i];
        }
    }
    return *this;
}
Filter::~Filter() {
    delete[] w;
}

Filter Filter::low(int n, int sr, double c) { 
    Filter f(n);
    for (int i = 0;i < n;i++) { 
        double freq = (i <= n / 2 ? i : n - i) * (double)sr / n;
        f.w[i] = (freq <= c ? 1 : 0);
    }
    return f;
}


Filter Filter::high(int n, int sr, double c) {
    Filter f(n);
    for (int i = 0;i < n;i++) {
        double freq = (i <= n / 2 ? i : n - i) * (double)sr / n;
        f.w[i] = (freq >= c ? 1 : 0); 
    }return f;
}
Filter Filter::band(int n, int sr, double l, double h) {
    Filter f(n);
    for (int i = 0;i < n;i++) {
        double freq = (i <= n / 2 ? i : n - i) * (double)sr / n;
        f.w[i] = (freq >= l && freq <= h ? 1 : 0); 
    }
    return f;
}

void bitReverse(Complex* a, int n) { 
    for (int i = 1, j = 0;i < n;i++) { 
        int bit = n >> 1;
        for (;j & bit;bit >>= 1)j ^= bit;j ^= bit;if (i < j) { 
            Complex t = a[i];
            a[i] = a[j];
            a[j] = t; 
        } 
    } 
}

void fft(Complex* a, int n, bool inv) {
    bitReverse(a, n);
    for (int len = 2;len <= n;len <<= 1)
    {
        double ang = 2 * PI / len * (inv ? -1 : 1);
        Complex wlen(cos(ang), sin(ang));
        for (int i = 0;i < n;i += len) {
            Complex w(1, 0);
            for (int j = 0;j < len / 2;j++) {
                Complex u = a[i + j];
                Complex v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v
                    ;w = w * wlen;
            }
        }
    }
    if (inv) {
        for (int i = 0;i < n;i++) {
            a[i].re /= n;a[i].im /= n;
        }
    }
}

Spectrum::Spectrum(int n) :n(n) {
    data = new Complex[n];
}
Spectrum::Spectrum(const Spectrum& other) : n(other.n), orig(other.orig), rate(other.rate), header(other.header) {
    data = new Complex[n];
    for (int i = 0;i < n;i++) {
        data[i] = other.data[i];
    }
}
Spectrum& Spectrum::operator=(const Spectrum& other) {
    if (this != &other) {
        delete[] data;
        n = other.n;
        orig = other.orig;
        rate = other.rate;
        header = other.header;
        data = new Complex[n];
        for (int i = 0;i < n;i++) {
            data[i] = other.data[i];
        }
    }
    return *this;
}
Spectrum::~Spectrum() {
    delete[] data;
}
void Spectrum::apply(Filter& f) {
    for (int i = 0;i < n;i++) {
        data[i].re *= f.w[i];
        data[i].im *= f.w[i]; 
    }
}

Whisper Spectrum::IFFT() {
    fft(data, n, true);
    Whisper w;
    w.header = header;
    w.n = orig;
    w.samples = new short[orig];
    for (int i = 0;i < orig;i++) {
        w.samples[i] = clip(data[i].re);
    }
    w.header.subchunk2Size = w.n * 2;
    w.header.chunkSize = 36 + w.header.subchunk2Size;
    return w;
}

Whisper::Whisper() :samples(NULL), n(0) {}

Whisper::Whisper(const char* path) {
    ifstream f(path, ios::binary);

    if (!f) {
        cout << "[ERROR] Cannot open file: " << path << endl;
        samples = NULL;
        n = 0;
        return;
    }

    f.read((char*)&header, sizeof(header));

    if (!f) {
        cout << "[ERROR] Invalid WAV file\n";
        samples = NULL;
        n = 0;
        return;
    }

    int bytes = header.subchunk2Size;
    int count = bytes / (header.bitsPerSample / 8);

    n = count;
    samples = new short[n];

    if (header.bitsPerSample == 8) {
        unsigned char x;
        for (int i = 0; i < n; i++) {
            f.read((char*)&x, 1);
            samples[i] = (x - 128) << 8;
        }
    }
    else {
        f.read((char*)samples, bytes);
    }

    header.bitsPerSample = 16;
    f.close();
}

Whisper::~Whisper() { 
    delete[] samples;
}

void Whisper::normalize() { 
    int m = 0;
    for (int i = 0;i < n;i++) {
        if (abs(samples[i]) > m)
            m = abs(samples[i]);
    }
    double s = 32767.0 / m;
    for (int i = 0;i < n;i++) {
        samples[i] = clip(samples[i] * s);
    }
}

Whisper& Whisper::operator*=(double g) { 
    for (int i = 0;i < n;i++) {
        samples[i] = clip(samples[i] * g);
    }
    return *this;
}

void Whisper::reverse() {
    for (int i = 0;i < n / 2;i++) {
        short t = samples[i];
        samples[i] = samples[n - i - 1];
        samples[n - i - 1] = t; 
    }
}

Whisper Whisper::operator+(const Whisper& o) { 
    Whisper w;
    w.n = (n < o.n ? n : o.n);
    w.samples = new short[w.n];
    for (int i = 0;i < w.n;i++) {
        w.samples[i] = clip(samples[i] + o.samples[i]);
    }
    w.header = header;
    w.header.subchunk2Size = w.n * 2;
    w.header.chunkSize = 36 + w.header.subchunk2Size;
    return w;
}

Whisper& Whisper::operator+=(const Whisper& o) {
    short int* ns = new short int[n + o.n];
    for (int i = 0;i < n;i++) {
        ns[i] = samples[i];
    }
    for (int i = 0;i < o.n;i++) {
        ns[n + i] = o.samples[i];
    }
    delete[] samples;
    samples = ns;
    n += o.n;
    header.subchunk2Size = n * 2;
    header.chunkSize = 36 + header.subchunk2Size;
    return *this;
}

Whisper Whisper::splice(double startSec, double endSec) {
    int startIdx = (int)(startSec * header.sampleRate);
    int endIdx = (int)(endSec * header.sampleRate);
    if (startIdx < 0) startIdx = 0;
    if (endIdx > n) endIdx = n;
    if (startIdx > endIdx) startIdx = endIdx;

    Whisper w;
    w.header = header;
    w.n = endIdx - startIdx;
    w.samples = new short[w.n];
    for (int i = 0;i < w.n;i++) {
        w.samples[i] = samples[startIdx + i];
    }
    w.header.subchunk2Size = w.n * 2;
    w.header.chunkSize = 36 + w.header.subchunk2Size;
    return w;
}


Spectrum Whisper::FFT() {
    int size = nextPow2(n);
    Spectrum s(size);
    s.orig = n;
    s.rate = header.sampleRate;
    s.header = header;
    for (int i = 0;i < size;i++) {
        s.data[i] = (i < n ? Complex(samples[i], 0) : Complex(0, 0));
    }
    fft(s.data, size, false);
    return s;
}

void Whisper::save(const char* path) { 
    header.subchunk2Size = n * 2;
    header.chunkSize = 36 + header.subchunk2Size;
    ofstream f(path, ios::binary);
    f.write((char*)&header, sizeof(header));
    f.write((char*)samples, n * 2);
    f.close(); 
}


void timeDomainMenu(Whisper* whispers[4], Spectrum* spectrums[4], int& activeW, int& activeS) {
    int choice;
    while (true) {
        cout << "\nTIME DOMAIN (WHISPER) MENU\n";
        cout << "1. Load Whisper\n";
        cout << "2. Set Active Whisper\n";
        cout << "3. Apply Gain\n";
        cout << "4. Normalize\n";
        cout << "5. Reverse\n";
        cout << "6. Splice\n";
        cout << "7. FFT -> Spectrum\n";
        cout << "8. Save Active Whisper\n";
        cout << "0. Back\n";
        cout << "Choice: ";
        cin >> choice;

        if (choice == 0) return;

        //launching whisper
        if (choice == 1) {
            int index;
            char path[100];

            cout << "Enter index (0-3): ";
            cin >> index;
            cout << "Enter file path: ";
            cin >> path;

            if (index >= 0 && index < 4) {
                delete whispers[index]; //avoid memory leak
                Whisper* temp = new Whisper(path);

                if (temp->samples == NULL) {
                    cout << "[ERROR] Loading failed. Slot unchanged.\n";
                    delete temp;
                }
                else {
                    delete whispers[index];
                    whispers[index] = temp;
                    cout << "[OK] Loaded into slot " << index << "\n";
                }
            }
            else {
                cout << "[ERROR] Invalid index\n";
            }
        }

        //bringing whisper active
        else if (choice == 2) {
            int index;
            cout << "Enter index (0-3): ";
            cin >> index;

            if (index >= 0 && index < 4 && whispers[index] != NULL) {
                activeW = index;
                cout << "[OK] Active Whisper = " << index << "\n";
            }
            else {
                cout << "[ERROR] Invalid or empty slot\n";
            }
        }

        //applying the gain
        else if (choice == 3) {
            if (activeW == -1) {
                cout << "[ERROR] No active Whisper\n";
                continue;
            }

            double g;
            cout << "Enter gain factor: ";
            cin >> g;

            (*whispers[activeW]) *= g;
            cout << "[OK] Gain applied\n";
        }

        //clamping between the range
        else if (choice == 4) {
            if (activeW == -1) {
                cout << "[ERROR] No active Whisper\n";
                continue;
            }

            whispers[activeW]->normalize();
            cout << "[OK] Normalized\n";
        }

        //bit reversal
        else if (choice == 5) {
            if (activeW == -1) {
                cout << "[ERROR] No active Whisper\n";
                continue;
            }

            whispers[activeW]->reverse();
            cout << "[OK] Reversed\n";
        }

        //splice
        else if (choice == 6) {
            if (activeW == -1) {
                cout << "[ERROR] No active Whisper\n";
                continue;
            }

            int index;
            double startSec, endSec;
            cout << "Store Whisper index (0-3): ";
            cin >> index;
            cout << "Enter start time (sec): ";
            cin >> startSec;
            cout << "Enter end time (sec): ";
            cin >> endSec;

            if (index >= 0 && index < 4) {
                delete whispers[index];
                whispers[index] = new Whisper(whispers[activeW]->splice(startSec, endSec));
                cout << "[OK] Spliced into Whisper slot " << index << "\n";
            }
            else {
                cout << "[ERROR] Invalid index\n";
            }
        }

        //fft
        else if (choice == 7) {
            if (activeW == -1) {
                cout << "[ERROR] No active Whisper\n";
                continue;
            }

            int index;
            cout << "Store Spectrum index (0-3): ";
            cin >> index;

            if (index >= 0 && index < 4) {
                delete spectrums[index];
                spectrums[index] = new Spectrum(whispers[activeW]->FFT());
                cout << "[OK] FFT stored in Spectrum slot " << index << "\n";
            }
            else {
                cout << "[ERROR] Invalid index\n";
            }
        }

        //save
        else if (choice == 8) {
            if (activeW == -1) {
                cout << "[ERROR] No active Whisper\n";
                continue;
            }

            char out[100];
            cout << "Enter output file path: ";
            cin >> out;

            whispers[activeW]->save(out);
            cout << "[OK] File saved\n";
        }

        else {
            cout << "[ERROR] Invalid choice\n";
        }
    }
}

void frequencyDomainMenu(Whisper* whispers[4], Spectrum* spectrums[4], int& activeW, int& activeS) {
    int choice;
    while (true) {
        cout << "\nFREQUENCY DOMAIN (SPECTRUM) MENU\n";
        cout << "1. Set Active Spectrum\n";
        cout << "2. Low Pass Filter\n";
        cout << "3. High Pass Filter\n";
        cout << "4. Band Pass Filter\n";
        cout << "5. IFFT -> Whisper\n";
        cout << "0. Back\n";
        cout << "Choice: ";
        cin >> choice;

        if (choice == 0) return;

        //activating spectrum
        if (choice == 1) {
            int index;
            cout << "Enter index (0-3): ";
            cin >> index;

            if (index >= 0 && index < 4 && spectrums[index] != NULL) {
                activeS = index;
                cout << "[OK] Active Spectrum = " << index << "\n";
            }
            else {
                cout << "[ERROR] Invalid or empty slot\n";
            }
        }

        //other factors filterations
        else if (choice == 2 || choice == 3 || choice == 4) {
            if (activeS == -1) {
                cout << "[ERROR] No active Spectrum\n";
                continue;
            }

            double c1, c2;

            if (choice == 2) {
                cout << "Enter cutoff frequency (Hz): ";
                cin >> c1;
                Filter f = Filter::low(spectrums[activeS]->n, spectrums[activeS]->rate, c1);
                spectrums[activeS]->apply(f);
            }

            else if (choice == 3) {
                cout << "Enter cutoff frequency (Hz): ";
                cin >> c1;
                Filter f = Filter::high(spectrums[activeS]->n, spectrums[activeS]->rate, c1);
                spectrums[activeS]->apply(f);
            }

            else {
                cout << "Enter low cutoff (Hz): ";
                cin >> c1;
                cout << "Enter high cutoff (Hz): ";
                cin >> c2;
                Filter f = Filter::band(spectrums[activeS]->n, spectrums[activeS]->rate, c1, c2);
                spectrums[activeS]->apply(f);
            }

            cout << "[OK] Filter applied\n";
        }

        //ifft
        else if (choice == 5) {
            if (activeS == -1) {
                cout << "[ERROR] No active Spectrum\n";
                continue;
            }

            int index;
            cout << "Store Whisper index (0-3): ";
            cin >> index;

            if (index >= 0 && index < 4) {
                delete whispers[index];
                whispers[index] = new Whisper(spectrums[activeS]->IFFT());
                cout << "[OK] IFFT stored in Whisper slot " << index << "\n";
            }
            else {
                cout << "[ERROR] Invalid index\n";
            }
        }

        else {
            cout << "[ERROR] Invalid choice\n";
        }
    }
}

int main() {
    //4 whisper and spectrum slots
    Whisper* whispers[4] = { NULL, NULL, NULL, NULL };
    Spectrum* spectrums[4] = { NULL, NULL, NULL, NULL };

    //active places
    int activeW = -1;
    int activeS = -1;

    int choice;

    while (true) {
        cout << "\nORCHESTRATOR SYSTEM\n";
        cout << "1. Time Domain (Whisper) Menu\n";
        cout << "2. Frequency Domain (Spectrum) Menu\n";
        cout << "0. Exit\n";
        cout << "Choice: ";
        cin >> choice;

        if (choice == 0) break;
        else if (choice == 1) timeDomainMenu(whispers, spectrums, activeW, activeS);
        else if (choice == 2) frequencyDomainMenu(whispers, spectrums, activeW, activeS);
        else cout << "[ERROR] Invalid choice\n";
    }

    //cleanup of memory
    for (int i = 0; i < 4; i++) {
        delete whispers[i];
        delete spectrums[i];
    }

    return 0;
}