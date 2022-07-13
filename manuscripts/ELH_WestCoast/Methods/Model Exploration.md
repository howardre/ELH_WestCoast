**Data**
[[Rockfish Recruitment & Ecosystem Assessment Survey]]

**Species**
[[Pacific sanddab]] 
- Access to both YOY and adult, only using YOY from RREAS
	- No length measurements before 2011

[[Pacific hake]]
- Access to only YOY from RREAS
	- Length measurements from 1994

**Models**
- Simple GAMs
	- Either used Tweedie or Gaussian
	- Either with or without variable coefficient terms
	- PDO, NPGO, or ONI
	- Formulations:
		lncatch ~ f(year) + s(lon, lat) + s(depth) + s(DOY)
		lncatch ~ f(year) + s(lon, lat) + s(depth) + s(DOY) + s(climate index)
		lncatch ~ f(year) + s(lon, lat) + s(depth) + s(DOY) + s(lat, lon, by = climate index)
	- Results
		- [[Pacific sanddab]]
			- Gaussian: Model with s(ONI), $R^2$ = 0.356, deviance = 37.5%
			- Tweedie: Model with s(lat, lon, by = NPGO), deviance = 59.9%
		- [[Pacific hake]]
			- Gaussian: Model with s(lat, lon, by = NPGO), $R^2$ = 0.232, deviance = 25.8%
			- Tweedie: Model with s(lat, lon, by = NPGO), deviance = 47.3%
- Length binned GAMs
	- Split the data into two length bins
		- Used histograms to look at number to split at
	- Compare split models to base model using all the data
		- Had to limit due to lack of length data in some years
	- Results
		- [[Pacific sanddab]]
			- Data from 2011 onward
		- [[Pacific hake]]
			- Data from 1994 onward
			- Base model with s(lat, lon, by = NPGO), deviance = 46.1%
			- 