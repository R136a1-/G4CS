"""
This program creates a fiducial line utilizing a sliding window of desired size and the IQR of contained in that window.

"your_file_here" must contain magnitudes in two different bands (to create a CMD), indexing will need to be adjusted

Includes optional return statement in the main function along with optional plots to display status of ridgeline

"""
import numpy as np
import collections as clt
from matplotlib import pyplot as plt
import math
from decimal import *
from scipy.stats import norm
from scipy.interpolate import interp1d
from matplotlib.ticker import MultipleLocator

class RidgeLine:
	def main(self):
		self.all_info = np.loadtxt("your_file_here")
		
		print "	"
		print "------------------CMD/RIDGELINE PARAMETERS-------------------"
		print "	"
		self.cmd_increments = int(raw_input("How many steps would you like to take? "))
		self.min_mag = float(raw_input("What is the faintest magnitude you would like to consider? "))
		self.max_mag = float(raw_input("What is the brightest magnitude you would like to consider? "))
		self.bin_size = Decimal(raw_input("What bin size would you like to have (in magnitudes) to smooth the ridgeline? "))

		self.ridgeline_colour, self.ridgeline_y, self.ridgeline_col_addsig, self.ridgeline_col_subsig, self.mag_filter = self.make_ridgeline()
		
		self.smooth_col, self.smooth_col_addsig, self.smooth_col_subsig, self.smooth_mags = self.smooth_ridgeline(self.bin_size)  #optional smoothing function to resample and smooth the fiducial line
		
		print "	"
		answer = raw_input("Would you like to try a different bin size and replot? (Y or N) ")

		while answer.lower() == "Y".lower():
			new_bin_size = Decimal(raw_input("Enter a new bin size: "))
			self.smooth_ridgeline(new_bin_size)
			answer = raw_input("Would you like to try a different bin size and replot? (Y or N) ")

		smooth_info = np.empty([self.smooth_mags.size, 2])

		smooth_info[:, 0] = self.smooth_col_addsig
		smooth_info[:, 1] = self.smooth_mags

		sort = np.argsort(smooth_info, 0)
	
		x = smooth_info[sort[:, 0], 0]
		y = smooth_info[sort[:, 1], 1]

		f = interp1d(x, y)

		return self.smooth_col, self.smooth_col_addsig, self.smooth_col_subsig, self.smooth_mags, f, self.min_mag, self.max_mag

	def make_ridgeline(self):
		"""
		This function first filters through the data keeping only data points within the desired magnitude range.
		It then creates "all_info_sorted" an array of col1 = F275W - Ks, col2 = F275W data points within mag range.

		"all_info_sorted" is then sorted by increasing magnitude in order to begin sliding the window.

		Iterating through "all_info_sorted", sliding a window of desired size through the array, a mean value of
		each group contained within the window is found, using 70% of the members.

		Return an array of the mean colour and mean magntitude contained within each sliding window.

		This is the ridgeline!
		""" 
		mag_filter = []

		for num in range(0, self.all_info.size/67):				#iterate through, only keep data within mag. window
			mag = self.all_info[num, 22]

			if self.max_mag < mag < self.min_mag:
				mag_filter.append(num)

		all_info_sorted = np.empty([len(mag_filter), 2]) 			#create new array to store colour and mag
		all_info_sorted[:, 0] = self.all_info[mag_filter, 22] - self.all_info[mag_filter, 11] 
		all_info_sorted[:, 1] = self.all_info[mag_filter, 22]
		all_info_sorted = all_info_sorted[all_info_sorted[:, 1].argsort()] 	#organize the new array in order of increasing Ks mag.

		good_colour = []
		good_colour_addsigma = []
		good_colour_subsigma = []
		good_band1_mag = []

		index1 = 0
		index2 = self.cmd_increments


		for num in range(0,	all_info_sorted.size/2):	 #iterate through, in desired window sizes calculate the mean 70% of the members in the group
			group =	all_info_sorted[index1:index2, :]

			sorted_group = group[group[:, 0].argsort()]
			sorted_group_70 = sorted_group[0:((group.size/2) * 0.70), :]

			try:
				mean_colour = sorted_group_70[int(sorted_group_70.size/4 - 1), 0]
				sigma_colour = 0.63 * math.fabs(sorted_group_70[int(3*sorted_group_70.size/8 - 1), 0] - sorted_group_70[int(sorted_group_70.size/8 - 1), 0]) #interquartial range, sigma may not be classical definition, may change with binary fraction
				mean_band1_mag = sorted_group_70[int(sorted_group_70.size/4 - 1), 1]									     #take 35th data point = median of first 69
																					     #take the difference between (17th dp and 52 dp)*0.63 = sigma
				good_colour.append(mean_colour)
				good_colour_addsigma.append(mean_colour + sigma_colour)
				good_colour_subsigma.append(mean_colour - sigma_colour)
				good_band1_mag.append(mean_band1_mag)

				index1 += 1
				index2 += 1

			except:
				break

	
		final_good_colour = np.delete(np.asarray(good_colour), len(good_colour) - 1, 0)
		final_good_col_addsig = np.delete(np.asarray(good_colour_addsigma), len(good_colour_addsigma) - 1, 0)
		final_good_col_subsig = np.delete(np.asarray(good_colour_subsigma), len(good_colour_subsigma) - 1, 0)
		final_good_band1_mag = np.delete(np.asarray(good_band1_mag), len(good_band1_mag) - 1, 0)

		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.scatter(self.all_info[mag_filter, 22] - self.all_info[mag_filter, 11], self.all_info[mag_filter, 22], marker = ".", color = "r", label="PM Cleaned Data")  
		ax1.scatter(final_good_colour, final_good_band1_mag, marker = ".", label="Filtered Data")
		ax1.scatter(final_good_col_addsig, final_good_band1_mag, marker = ".", color= "y", label="Filtered Data + Sigma")
		ax1.scatter(final_good_col_subsig, final_good_band1_mag, marker = ".", color="b", label="Filtered Data - Sigma")

		plt.xlim(0, 8)
		plt.ylim(24, 14)
		plt.show()

		return final_good_colour, final_good_band1_mag, final_good_col_addsig, final_good_col_subsig, mag_filter


	def smooth_ridgeline(self, bin_size):
		binned_info = clt.defaultdict(list)
		binned_addsig = clt.defaultdict(list)
		binned_subsig = clt.defaultdict(list)
		binned_means = {}

		bin_sig_figs = int(str(bin_size.as_tuple().exponent).strip("-"))
		getcontext().prec = bin_sig_figs
		sf = bin_sig_figs + 3

		bin_start = float("{:.{}}".format(str(np.amin(self.ridgeline_y)), sf))
		bin_end = float("{:.{}}".format(str(np.amax(self.ridgeline_y)), sf))
		bin_size = float(bin_size)

		bins = np.arange(bin_start, bin_end, bin_size)

		for num in range(self.ridgeline_y.size):
			ridge_mag = self.ridgeline_y[num]
			ridge_colour = self.ridgeline_colour[num]
			ridge_col_addsig = self.ridgeline_col_addsig[num]
			ridge_col_subsig = self.ridgeline_col_subsig[num]

			for num2 in range(bins.size):
				bin_mag = bins[num2]

				if math.fabs(ridge_mag - bin_mag) <= bin_size:
					binned_info[bin_mag].append(ridge_colour)
					binned_addsig[bin_mag].append(ridge_col_addsig)
					binned_subsig[bin_mag].append(ridge_col_subsig)


		mean_mags = []
		mean_col = []
		mean_col_addsig = []
		mean_col_subsig = []
		sigma_col = []

		for num in range(len(binned_info.keys())):
			col = np.mean(np.asarray(binned_info[binned_info.keys()[num]]))
			col_addsig = np.mean(np.asarray(binned_addsig[binned_addsig.keys()[num]]))
			col_subsig = np.mean(np.asarray(binned_subsig[binned_subsig.keys()[num]]))

			sig = np.std(np.asarray(binned_info[binned_info.keys()[num]]))

			binned_means[binned_info.keys()[num]] = col
			mean_mags.append(binned_info.keys()[num])
			mean_col.append(col)
			mean_col_addsig.append(col_addsig)
			mean_col_subsig.append(col_subsig)
			sigma_col.append(sig)


		mean_mags = np.asarray(mean_mags)
		mean_col = np.asarray(mean_col)
		mean_col_addsig = np.asarray(mean_col_addsig)
		mean_col_subsig = np.asarray(mean_col_subsig)
		sigma_col = np.asarray(sigma_col)

		fig = plt.figure()
		ax1 = fig.add_subplot(111)
		ax1.scatter(self.all_info[self.mag_filter, 22] - self.all_info[self.mag_filter, 11], self.all_info[self.mag_filter, 22], marker = ".", color = "r", label="PM Cleaned Data")  
		#ax1.errorbar(mean_col, mean_mags, xerr=sigma_col, fmt=".", color = "k", label="Smoothed Ridgeline")
		#ax1.scatter(mean_col, mean_mags, marker=".")
		ax1.scatter(mean_col_addsig, mean_mags, marker=".")
		ax1.scatter(mean_col_subsig, mean_mags, marker=".")

		
		plt.xlim(-2, 8)
		plt.ylim(self.min_mag + 2, self.max_mag)
		plt.legend(loc='upper left')
		plt.xlabel(r"$F336W-K_{s}$")
		plt.ylabel(r"$F336W$")
		plt.title("Smoothed Ridgeline")
		plt.show()

		return mean_col, mean_col_addsig, mean_col_subsig, mean_mags


if __name__ == "__main__":
    X = RidgeLine()
    X.main()
