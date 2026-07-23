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



* **s03\_participants\_demographics:** Jasp-file to calculate descriptive information, contingency tables with chi-squared tests and Mann-Whitney U tests for group differences



* **s04\_partial\_correlations\_epoc:** Jasp-file to calculate partial correlation analyses for each questionnaire score while controlling for the respective other questionnaires, age, and education and save it in EPOC\\03\_figures\\06\_heatmaps



* **s05\_partial\_correlations\_covidom:** Jasp-file calculating the same analyses from s05, but for the covidom dataset



* **s06\_pre\_post\_increase\_eeg:** R-script to look for outliers in the EEG data and calculate graphics for pre-post changes in alpha and theta Peak power (graphic saved in EPOC\\03\_figures\\06\_eeg\_peaks)

&#x20;

* **s07\_EEG\_correlations:** Jasp-file to calculate Wilcoxon signed-rank tests for differences pre- vs post-stimulus, and Mann-Whitney U tests for group differences in the increase of EEG peak power between particpants with clinically significant fatigue and participants without clinically significant fatigue, for alpha and Theta
* 
* **s08\_topoplots\_beforeafter\_groups**: Calculating Time-Frequency Analysis on data and creating topographical plots for Alpha and Theta seperated by Group and pre- vs post-stimulus (in EPOC\\03\_figures\\06\_eeg\_peaks) + saving data in EPOC\\02\_data\\01\_prep\\aperiodic\_pvt\\topoplot\_workspace



* **s09\_TFR\_plots**: creating TFR plots based on data from s08\_topoplots\_beforeafter\_groups in EPOC\\03\_figures\\06\_eeg\_peaks



* **s10\_subjective\_cognitive\_symptoms:** Jasp-file to calculate Mann-Whitney U tests for group differences in PVT reaction time between participants reporting cognitive symptoms and participants without such symptoms



* **s11\_mediation\_epoc:** R-script to calculate a regression-based mediation analysis, testing the indirect effect of fatigue in the association between subjective cognitive symptoms and PVT reaction time, and the remaining direct effect
* 
* **s12\_holm\_corrections:** R-script to calculate holm-corrected p-values 

