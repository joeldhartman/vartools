library(tseries)
library(forecast)

DoArimaFitPlot <- function(mag, plotoutdir, lcbasename) {
	mag_ts <- ts(mag, start=1, end=length(mag), frequency=1)
	arima_model <- auto.arima(mag_ts)
	mag_arima <- mag - as.vector(arima_model$residuals)
	mag_forecast <- forecast(arima_model)
	png(paste(plotoutdir, lcbasename, ".arimaforecast.png", sep=""), width = 640, height=480)
	plot(mag_forecast)
	dev.off()
	png(paste(plotoutdir, lcbasename, ".arimaresiduals.png", sep=""), width = 640, height=480)
	checkresiduals(arima_model)
	dev.off()
	return(mag_arima)
}
