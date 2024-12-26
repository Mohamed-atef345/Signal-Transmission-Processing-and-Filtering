# Communication System

This project simulates a communication system using MATLAB, where you can process audio signals, simulate channel effects, add noise, and filter the signal at the receiver. It provides a hands-on demonstration of how signals are transmitted and processed in real-world scenarios.

## Features
- Visualizes the audio signal in time and frequency domains.
- Simulates channel effects with various predefined channels.
- Adds Gaussian noise to the signal with user-defined variance.
- Filters the noisy signal to retrieve the original signal.
- Supports custom audio files.

## Prerequisites
- MATLAB installed on your system.

## How to Use
1. Clone this repository to your local machine.
   ```bash
   git clone <repository-url>
   ```
2. Open MATLAB and navigate to the project directory.
3. Add your audio file in the same directory as the code. Ensure the file is in `.mp3` format.
4. Update the filename in the code to match your audio file. By default, the code is set to `audio_file.mp3`. For example:
   ```matlab
   [x, fs] = audioread('your_audio_file_name.mp3');
   ```
5. Run the `communication_system.m` file in MATLAB.
6. Follow the on-screen instructions to interact with the program.

## Code Workflow
1. **Transmission**:
   - Reads and plays the specified audio file.
   - Visualizes the signal in time and frequency domains.
2. **Channel Effects**:
   - Simulates the signal passing through predefined channels.
   - Visualizes the altered signals.
3. **Noise Addition**:
   - Adds Gaussian noise with user-specified variance.
   - Visualizes noisy signals.
4. **Receiver**:
   - Filters unwanted frequencies.
   - Converts the signal back to the time domain.
   - Visualizes and plays the filtered signals.

