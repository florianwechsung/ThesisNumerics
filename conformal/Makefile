.PHONY: levelset, levelset_plot, levelset_pv, stokes, stokes_pv, stokes_cdf, annulus, annulus_paraview

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
	${PVPYTHON} levelset_paraview.py

stokes:
	python3 stokes_conformal.py --base_inner elasticity --use_cr 1 
	python3 stokes_conformal.py --base_inner elasticity --use_cr 0 
	python3 stokes_conformal.py --base_inner laplace    --use_cr 1 
	python3 stokes_conformal.py --base_inner laplace    --use_cr 0 

stokes_pv:
	${PVPYTHON} stokes_paraview.py --base_inner elasticity --use_cr 1 
	${PVPYTHON} stokes_paraview.py --base_inner elasticity --use_cr 0 
	${PVPYTHON} stokes_paraview.py --base_inner laplace    --use_cr 1 
	${PVPYTHON} stokes_paraview.py --base_inner laplace    --use_cr 0 

stokes_cdf:
	python3 stokes_inv_cdf.py

annulus:
	#gmsh -2 -clscale 0.6 -smooth 100 annulus.geo
	python3 annulus_conformal.py --maxiter 50 --base laplace --rstar 0.5
	python3 annulus_conformal.py --maxiter 50 --base laplace --rstar 0.7
	python3 annulus_conformal.py --maxiter 50 --base laplace --rstar 0.8
	python3 annulus_conformal.py --maxiter 50 --base laplace --alpha 0.01 --rstar 0.7
	python3 annulus_conformal.py --maxiter 50 --base laplace --alpha 0.01 --rstar 0.79
	python3 annulus_conformal.py --maxiter 50 --base elasticity --rstar 0.79

annulus_pv:
	${PVPYTHON} annulus_paraview.py --base laplace --rstar 0.5
	${PVPYTHON} annulus_paraview.py --base laplace --rstar 0.7
	${PVPYTHON} annulus_paraview.py --base laplace --rstar 0.8
	${PVPYTHON} annulus_paraview.py --base laplace --alpha 0.01 --rstar 0.7
	${PVPYTHON} annulus_paraview.py --base laplace --alpha 0.01 --rstar 0.79
	${PVPYTHON} annulus_paraview.py --base elasticity --rstar 0.79
