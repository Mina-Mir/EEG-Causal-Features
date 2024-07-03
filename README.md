## Neural Oscillation Discovery for Working Memory Performance

[Fig.1.pdf](https://github.com/user-attachments/files/16090955/Fig.1.pdf)

This repository contains MATLAB functions designed for discovering neural oscillations that show causal relationships with working memory performance using EEG data.

In addition to my own implementations, I have integrated code from the following sources:
- [SingleCI Causal Feature Selection](https://gitlab.tuebingen.mpg.de/amastakouri/singleCICausalFeatureSelection)
- [Causal Learner Toolbox](https://github.com/z-dragonl/Causal-Learner)

### Usage

extract_performance.m function calculates the performance during N-back task using EEG data
power_calculation.m function calculates the power of neural oscillations from pre- and post-stimulus time windows
causal_features.m calculates the causal relations between the performance and power from neural oscillations

To use these functions, clone the repository and follow the instructions in each module's README file for setup and execution.

### References

If you use these methods in your research, please cite the respective sources appropriately.

