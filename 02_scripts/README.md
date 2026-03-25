## script description



* **s01\_badChannelsAndICA**: Read and filter data specifically for ICA Analysis; saving the following Information:

 	\* how many and which channels, segments and ICA components are considered noisy per participant in EPOC\\01\_data\\01\_individual\_data\_after\_ICA and Overall in EPOC\\01\_data\\noisy\_objects\_after\_s01.mat

 	\* figures showing z-scores of all segments in a scatterplot for each participant in EPOC\\03\_figures\\01\_bad\_trials\_z\_pvt

 	\* figures showing all components which are considered noisy and all components which will be kept in EPOC\\03\_figures\\02\_ica\_plots\_pvt



* **s02\_aperiodic\_power**: Using idividual files with noisy channels and components from EPOC\\01\_data\\individual\_data\_after\_ICA to calculate powerspectra and seperate the aperiodic component, resulting in

 	\* figures showing the original powerspectrum and its aperiodic component (left side) and single Trial powerpectra + average of every participant for pre- and post-stimulus in EPOC\\03\_figures\\03\_fooof\_pvt

 	\* figures showing Alpha and Theta Peaks of the amplitudes for pre- and post-stimulus in EPOC\\03\_figures\\04\_alpha\_peak\_pvt and EPOC\\03\_figures\\05\_theta\_peak\_pvt 

 	\* saving peak amplitude and its frequency + mean amplitude for alpha and theta pre- and post-Stimulus in EPOC\\01\_data\\02\_fooof\_pvt\\fooofed.mat

 	\* saving trials which were rejected by pop\_jointprob in EPOC\\01\_data\\bad\_trials\_fooof\_s02



* **s03\_participants\_demographics:** Jasp-file to calculate descriptive information, contingency tables with chi-squared Tests and Mann-Whitney U tests for group differences 



* **s04\_heatmap\_behavioral:** R-Script to calculate heatmap for questionnaire scores, cognitive test scores, age, and education and save it in EPOC\\03\_figures\\06\_heatmaps



* **s05\_SEM\_epoc:** R-Script to calculate structural equation model with EPOC data



* **s06\_PVT\_group:** Jasp-file to calculate Mann-Whitney U test for group differences in PVT reaction time



* **s07\_SEM\_covidom:** R-Script to calculate structural equation model with COVIDOM data



* **s08\_heatmap\_scatterplots\_eeg:** R-Script to calculate heatmap for questionnaire scores, cognitive test scores, and EEG peak power (in EPOC\\03\_figures\\06\_heatmaps), scatterplots for the significant correlations (in EPOC\\03\_figures\\07\_eeg\_peaks), and graphics for pre-post changes in alpha and theta Peak power



* **s09\_participants\_eeg\_peaks:** Jasp-file to calculate Wilcoxon signed-rank tests for differences pre- vs post-stimulus, and Mann-Whitney U tests for group differences in EEG peak power, for alpha and theta 



* **s10\_topoplots\_beforeafter\_groups**: Calculating Time-Frequency Analysis on data and creating topographical plots for Alpha and Theta seperated by Group and pre- vs post-stimulus (in EPOC\\03\_figures\\07\_eeg\_peaks) + saving data in EPOC\\02\_data\\01\_prep\\aperiodic\_pvt\\topoplot\_workspace



* **s11\_TFR\_plots**: creating TFR plots based on data from s10\_topoplots\_beforeafter\_groups in EPOC\\03\_figures\\07\_eeg\_peaks
