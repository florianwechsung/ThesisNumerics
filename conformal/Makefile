.PHONY: levelset, levelset_plot, levelset_pv

#PVPYTHON=/home/wechsung/bin/ParaView-5.6.0-MPI-Linux-64bit/bin/pvpython
PVPYTHON=/home/wechsung/bin/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython
	
levelset:
	python3 levelset_conformal.py --base-inner elasticity --clscale .1 --maxiter 100
	python3 levelset_conformal.py --base-inner elasticity --clscale .1 --maxiter 100 --alpha 0.3
	python3 levelset_conformal.py --base-inner elasticity --clscale .1 --maxiter 100 --alpha 0.01
	python3 levelset_conformal.py --base-inner elasticity --clscale .1 --maxiter 100 --alpha 0.001
	python3 levelset_conformal.py --base-inner laplace --clscale .1 --maxiter 100
	python3 levelset_conformal.py --base-inner laplace --clscale .1 --maxiter 100 --alpha 0.3
	python3 levelset_conformal.py --base-inner laplace --clscale .1 --maxiter 100 --alpha 0.01
	python3 levelset_conformal.py --base-inner laplace --clscale .1 --maxiter 100 --alpha 0.001

levelset_plot:
	python3 levelset_plot.py

levelset_pv:
	${PVPYTHON} levelset_pv.py
