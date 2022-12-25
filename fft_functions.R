library(ggplot2)
library(dplyr)

waveform <- function(freq_hzs, amplitudes, decay_tau = NA, duration_s = 1.0, sr = 200) {
  t_samples <- seq(sr * duration_s)
  waveform_matrix = matrix(nrow = length(freq_hzs), ncol = length(t_samples))
  
  for (i in seq_along(freq_hzs)) {
    if (is.na(decay_tau)) {
      waveform_matrix[i, ] <- amplitudes[i] * cos(2 * pi * freq_hzs[i] * t_samples / sr)
    }
    else {
      waveform_matrix[i, ] <- amplitudes[i] * cos(2 * pi * freq_hzs[i] * t_samples / sr) * exp(-t_samples / (decay_tau * sr))
    }
  }
  
  wv <- colSums(waveform_matrix)
  
  data.frame(
    t_samples = t_samples,
    waveform = wv,
    waveform_norm = wv / length(wv)
  )
}

ff_transform <- function(waveform_df) {
  original_fft <- fft(waveform_df$waveform_norm)
  y <- original_fft[1:floor(length(original_fft) / 2)]
  
  data.frame(
    fft_idx = seq_along(y),
    fft_complex = y,
    fft_re = Re(y),
    fft_im = Im(y),
    fft_mod = Mod(y)
  )
}

waveform_ggplot <- function(waveform_df) {
  ggplot(waveform_df, aes(x = t_samples, y = waveform_norm)) +
    geom_line() +
    xlim(c(0, nrow(waveform_df)))
}

fft_ggplot <- function(fft_df, fft_mod_cutoff = 1e-2) {
  fft_df %>%
    filter(fft_mod > fft_mod_cutoff) %>%
    ggplot() +
    geom_segment(aes(x = fft_idx, y = 0, xend = fft_idx, yend = fft_mod)) +
    xlim(c(0, nrow(fft_df)))
}

