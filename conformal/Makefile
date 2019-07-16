.PHONY: levelset, levelset_plot, levelset_pv, levelset_cdf, stokes, stokes_pv, stokes_cdf, annulus, annulus_paraview

PVPYTHON=/home/wechsung/bin/ParaView-5.6.0-osmesa-MPI-Linux-64bit/bin/pvpython
	
levelset:
	python3 levelset_conformal.py --base-inner elasticity --clscale .1 --maxiter 100
	python3 levelset_conformal.py --base-inner elasticity --clscale .1 --maxiter 100 --alpha 0.3
	python3 levelset_conformal.py --base-inner elasticity --clscale .1 --maxiter 100 --alpha 0.01
	python3 levelset_conformal.py --base-inner elasticity --clscale .1 --maxiter 100 --alpha 0.001
	python3 levelset_conformal.py --base-inner laplace --clscale .1 --maxiter 30
	python3 levelset_conformal.py --base-inner laplace --clscale .1 --maxiter 100 --alpha 0.3
	python3 levelset_conformal.py --base-inner laplace --clscale .1 --maxiter 100 --alpha 0.01
	python3 levelset_conformal.py --base-inner laplace --clscale .1 --maxiter 100 --alpha 0.001

levelset_convergence:
	python3 levelset_convergence.py

levelset_cdf:
	python3 levelset_cdf.py

levelset_pv:
	${PVPYTHON} levelset_paraview.py

stokes:
	gmsh -2 -clscale 0.2 -smooth 100 Sphere2D.geo
	python3 stokes_conformal.py --base_inner elasticity --use_cr 1 &
	python3 stokes_conformal.py --base_inner elasticity --use_cr 0 &
	python3 stokes_conformal.py --base_inner laplace    --use_cr 1 &
	python3 stokes_conformal.py --base_inner laplace    --use_cr 0 &

stokes_pv:
	${PVPYTHON} stokes_paraview.py --base_inner elasticity --use_cr 1 
	${PVPYTHON} stokes_paraview.py --base_inner elasticity --use_cr 0 
	${PVPYTHON} stokes_paraview.py --base_inner laplace    --use_cr 1 
	${PVPYTHON} stokes_paraview.py --base_inner laplace    --use_cr 0 

stokes_cdf:
	python3 stokes_cdf.py

stokes_convergence:
	python3 stokes_convergence.py

annulus:
	gmsh -2 -clscale 1.0 -smooth 100 annulus.geo
	python3 annulus_conformal.py --maxiter 50 --base laplace --rstar 0.75
	python3 annulus_conformal.py --maxiter 50 --base laplace --rstar 0.75 --alpha 0.01
	python3 annulus_conformal.py --maxiter 50 --base elasticity --rstar 0.75
	python3 annulus_conformal.py --maxiter 50 --base elasticity --rstar 0.75 --alpha 0.01
	python3 annulus_conformal.py --maxiter 50 --base laplace --rstar 0.80
	python3 annulus_conformal.py --maxiter 50 --base laplace --rstar 0.80 --alpha 0.01
	python3 annulus_conformal.py --maxiter 50 --base elasticity --rstar 0.80
	python3 annulus_conformal.py --maxiter 50 --base elasticity --rstar 0.80 --alpha 0.01
	python3 annulus_conformal.py --maxiter 50 --base laplace --rstar 0.85
	python3 annulus_conformal.py --maxiter 50 --base laplace --rstar 0.85 --alpha 0.01
	python3 annulus_conformal.py --maxiter 50 --base elasticity --rstar 0.85
	python3 annulus_conformal.py --maxiter 50 --base elasticity --rstar 0.85 --alpha 0.01

annulus_pv:
	${PVPYTHON} annulus_paraview.py --base laplace --rstar 0.75
	${PVPYTHON} annulus_paraview.py --base laplace --rstar 0.75 --alpha 0.01
	${PVPYTHON} annulus_paraview.py --base elasticity --rstar 0.75
	${PVPYTHON} annulus_paraview.py --base elasticity --rstar 0.75 --alpha 0.01
	${PVPYTHON} annulus_paraview.py --base laplace --rstar 0.80
	${PVPYTHON} annulus_paraview.py --base laplace --rstar 0.80 --alpha 0.01
	${PVPYTHON} annulus_paraview.py --base elasticity --rstar 0.80
	${PVPYTHON} annulus_paraview.py --base elasticity --rstar 0.80 --alpha 0.01
	${PVPYTHON} annulus_paraview.py --base laplace --rstar 0.85
	${PVPYTHON} annulus_paraview.py --base laplace --rstar 0.85 --alpha 0.01
	${PVPYTHON} annulus_paraview.py --base elasticity --rstar 0.85
	${PVPYTHON} annulus_paraview.py --base elasticity --rstar 0.85 --alpha 0.01
