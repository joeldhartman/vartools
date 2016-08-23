def plotlc(lcname,outdir,t,ph,mag,P):
	import matplotlib.pyplot as plt
	lcbasename = lcname.split('/')[-1]
	plt.ioff()
	plt.figure(1)
	plt.subplot(211)
	plt.gca().invert_yaxis()
	tcorr = t - t[0]
	plt.plot(tcorr, mag, 'bo', markersize=0.5)
	plt.ylabel('magnitude')
	plt.title(lcbasename+' P='+str(P))
	plt.xlabel('time - '+str(t[0]))
	plt.subplot(212)
	plt.gca().invert_yaxis()
	plt.plot(ph, mag, 'bo', markersize=0.5)
	plt.ylabel('magnitude')
	plt.xlabel('phase')
	plt.savefig(outdir+'/'+lcbasename+'.png',format="png")
	plt.close(1)
