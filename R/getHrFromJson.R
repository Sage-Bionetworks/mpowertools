####### MAIN
#' extracts hr for each color band from avg pixel value per frame of video (processed hr) JSON data file
#'
#'
#' @param hrJsonFileLoc path to hr json file
#' @return list containing hr and confidence of the estimate for each color (red, green, blue)
#' @export
#' @examples
#' @author Meghasyam Tummalacherla 


#############################################################
# Wrapper function to take in json and give HR per color channel
#############################################################

getHrFromJson <- function(hrJsonFileLoc=NA, windowLen = 10, freqRange = c(0.5,4), bpforder = 128){
  
  #############################################################
  # Main Code Block
  #############################################################
  
  # If no json file exists
  dat1 = list(red = NA, green = NA, blue = NA, error = NA, samplingRate = NA)
  if(is.na(hrJsonFileLoc)){dat1$error = 'No JSON file'; return(dat1) }
  
  # Get HR data
  dat = tryCatch({ jsonlite::fromJSON(as.character(hrJsonFileLoc)) }, 
                 error = function(e){ NA })
  if(all(is.na(dat))){dat1$error = 'JSON file read error'; return(dat1) }
  
  # Get sampling rate
  samplingRate = tryCatch({ length(dat$timestamp)/(dat$timestamp[length(dat$timestamp)] - dat$timestamp[1]) }, 
                          error = function(e){ NA })
  if(!is.finite(samplingRate)){dat1$error = 'Sampling Rate calculated from timestamp is Inf or NaN / timestamp not found in json'; return(dat1) }
  if(samplingRate < 55){if(samplingRate > 22){bpforder = 64}else{bpforder = 32}}
  dat1$samplingRate = samplingRate
  # Convert window length from seconds to samples
  windowLen = round(samplingRate*windowLen)
  
  # Apply pre processing filter signal between freqRange
  mforder = 4*round(60*samplingRate/220) + 1 # order for the running mean based filter
  
  # Split each color into segments based on windowLen
  dat = tryCatch({ dat %>% dplyr::select(red, green, blue) %>% na.omit() %>% 
      lapply(window_signal, windowLen, ovlp=0.5) }, 
      error = function(e){ NA })
  if(all(is.na(dat))){dat1$error = 'red, green, blue cannot be read from JSON'; return(dat1) }
  
  # Get HR for each filtered segment of each color
  dat <- dat %>% lapply(function(dfl){
    dfl = tryCatch({
      apply(dfl,2,getHR,samplingRate)}, error = function(e){ NA })
    dfl = as.data.frame(t(dfl))
    colnames(dfl) = c('hr','confidence')
    return(dfl)
  })
  if(all(is.na(dat))){dat1$error = 'HR calculation error'; return(dat1) }
  
  dat$error = 'none'
  if(samplingRate < 55){dat$error = 'Low Sampling Rate,atleast 55FPS needed'}
  dat$samplingRate = samplingRate
  return(dat)
  
}

#############################################################
# Required Sub Functions
#############################################################

# Compute start/end timestamps for each window
window_start_end_times <- function(t, window_length, overlap) {
  seq_length <- length(t)
  if (seq_length < window_length) {
    window_length <- seq_length
    overlap <- 1
  }
  start_indices <- seq(1, seq_length, window_length * overlap)
  end_indices <- seq(window_length, seq_length, window_length * overlap)
  start_indices <- start_indices[1:length(end_indices)]
  start_times <- t[start_indices]
  end_times <- t[end_indices]
  window_start_end_times <- dplyr::tibble(
    window = seq(1, length(start_indices)),
    window_start_time = start_times,
    window_end_time = end_times,
    window_start_index = start_indices,
    window_end_index = end_indices)
  return(window_start_end_times)
}

# Split given signal into windows of given length
window_signal <- function(values, window_length = 256, overlap = 0.5) {
  start_end_times <- window_start_end_times(
    values, window_length = window_length, overlap = overlap)
  nstart <- start_end_times$window_start_index
  nend <- start_end_times$window_end_index
  wn <- rep(1, window_length)
  a <- apply(cbind(nstart, nend), 1, function(x, a, wn) {
    a[seq(x[1], x[2], 1)] * wn
  }, values, wn)
  colnames(a) <- 1:dim(a)[2]
  return(a)
}


# An optimizing wrapper around getHrFromTimeSeries
getHR <- function(x, samplingRate, minHR = 40, maxHR = 200){
  # Apply pre processing filter signal between freqRange
  mforder = 4*round(60*samplingRate/220) + 1 # order for the running mean based filter
  
  y <- getfilteredsignal(x, mforder = mforder,samplingRate = samplingRate,
                         freqRange = c(0.5,4))
  yHR <- getHrFromTimeSeries(y, samplingRate = samplingRate,
                             minHR = minHR, maxHR = maxHR)
  
  if(yHR[2] < 0.5){
    y <- getfilteredsignal(x, mforder = mforder,
                           freqRange = c(0.5,2), samplingRate = samplingRate)
    yHR <- getHrFromTimeSeries(y, samplingRate = samplingRate,
                               minHR = 40, maxHR = 120)
    if(yHR[2] < 0.5){
      y <- getfilteredsignal(x, mforder = mforder,
                             freqRange = c(1,4), samplingRate = samplingRate)
      yHR <- getHrFromTimeSeries(y, samplingRate = samplingRate,
                                 minHR = 60, maxHR = 240)
    }
  }
  
  if(yHR[2] < 0.5){
    mforder = 2*round(60*samplingRate/220) + 1 # order for the running mean based filter
    
    y <- getfilteredsignal(x, mforder = mforder,samplingRate = samplingRate,
                           freqRange = c(0.5,4))
    yHR <- getHrFromTimeSeries(y, samplingRate = samplingRate,
                               minHR = minHR, maxHR = maxHR)
    
    if(yHR[2] < 0.5){
      y <- getfilteredsignal(x, mforder = mforder,
                             freqRange = c(0.5,2), samplingRate = samplingRate)
      yHR <- getHrFromTimeSeries(y, samplingRate = samplingRate,
                                 minHR = 40, maxHR = 120)
      if(yHR[2] < 0.5){
        y <- getfilteredsignal(x, mforder = mforder,
                               freqRange = c(1,4), samplingRate = samplingRate)
        yHR <- getHrFromTimeSeries(y, samplingRate = samplingRate,
                                   minHR = 60, maxHR = 240)
      }
    }
  }
  
  return(yHR)
}

# Given a time series, get HR
getHrFromTimeSeries <- function(x, samplingRate, minHR = 40, maxHR=200){
  x[is.na(x)] <- 0
  
  
  spectX <- get_spectrum(x, sampling_rate = samplingRate, nfreq = 512) %>% 
    dplyr::filter(freq >= 0.5)
  estHR <- spectX$freq[which.max(spectX$pdf)]*samplingRate
  minHR <- max(estHR-12, minHR)
  maxHR <- min(estHR+12, maxHR)
  
  x <- stats::acf(x,lag.max = 250, plot=F)$acf
  y <- 0*x
  
  
  
  y[round(60*samplingRate/maxHR):round(60*samplingRate/minHR)] = x[round(60*samplingRate/maxHR):round(60*samplingRate/minHR)]
  confidence = max(y)/max(x)
  hr = 60*samplingRate/(which.max(y)-1)
  
  # If hr or condidence is NaN, then return hr = 0 and confidence = 0
  if(is.na(confidence) || is.na(hr)){
    confidence = NA
    hr = NA
  }
  
  return(c(hr, confidence))
}

# Bandpass and sorted mean filter the given signal
getfilteredsignal <- function(x, mforder = 65, bpforder = 128, freqRange=c(0.4,4), samplingRate){
  
  # Defaults are set for 60Hz sampling rate
  x[is.na(x)]<-0
  x = x-mean(x) #Mean centering the signal
  
  # Bandpass filter the given time series data
  if(samplingRate > 2*freqRange[2]){ 
    bandPassFilt = signal::fir1(bpforder-1, c(freqRange[1] * 2/samplingRate, freqRange[2] * 2/samplingRate),
                                type="pass", 
                                window = seewave::hamming.w(bpforder))
  }else{
    bandPassFilt = signal::fir1(bpforder, freqRange[1] * 2/samplingRate, # the order for a high pass filter needs to be even
                                type="high", 
                                window = seewave::hamming.w(bpforder+1))   
  }
  
  x = signal::filtfilt(bandPassFilt, x)
  x = x[((bpforder/2)+1):(length(x)-(bpforder/2)+1)]
  
  # Sorted Mean filter the given signal
  y = 0*x
  for (i in seq((mforder+1)/2, length(x)-(mforder-1)/2,1)){
    tempseq <- (x[(i-((mforder-1)/2)):((i+((mforder-1)/2)))])
    # y[i] = x[i]-(sum(tempseq)-max(tempseq))/(mforder-1)
    
    tempseq <- tempseq - mean(tempseq)
    # y[i] = (((x[i] - max(tempseq) - min(tempseq))-(sum(tempseq)-max(tempseq))/(mforder-1))/(max(tempseq)-min(tempseq) + 0.0000001))
    y[i] = (((x[i] - max(tempseq) - min(tempseq))-(sum(tempseq)-max(tempseq))/(mforder-1))/(max(tempseq)-min(tempseq) + 0.0000001))
  }
  y = y[((mforder+1)/2): (length(x)-(mforder-1)/2)]
  return(y)
}


get_spectrum <- function(values, sampling_rate = 60, nfreq = 512){
  tmp <- stats::spec.ar(values, n.freq = nfreq, plot = F)
  spectrum <- data.frame(freq = tmp$freq * sampling_rate, pdf = tmp$spec)
  return(spectrum)
}
