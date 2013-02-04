#a wrapper to make singlet follow the gating method convention
singlet.gate <- function(fs
						, xChannel = "FSC-A"
						,yChannel = "FSC-H"
						, prediction_level = 0.99
						,...) {
  
    # Creates a list of polygon gates based on the prediction bands at the minimum and maximum
    # x_channel observation using a robust linear model trained by flowStats.

	
        singletGate(fs[[1]], area = xChannel
					, height = yChannel
					,prediction_level = prediction_level
					)$gate
  
}

