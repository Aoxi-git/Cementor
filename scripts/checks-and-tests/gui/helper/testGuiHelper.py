#!/usr/bin/python
# -*- coding: utf-8 -*-
##########################################################################
#  2019        Janek Kozicki                                             #
#                                                                        #
#  This program is free software; it is licensed under the terms of the  #
#  GNU General Public License v2 or later. See file LICENSE for details. #
##########################################################################

import subprocess, time, os, sys
import yade
from yade import qt


class TestGUIHelper:
	"""
	This simple class makes screenshots.
	"""

	def __init__(self, name=None, recenter=False, cameraParams=None):
		self.recenter = recenter  # recenter after adjusting camera
		# cameraParams allows to use custom values of (v.lookAt, v.viewDir, v.eyePosition, v.upVector) # ←← in that order.
		# Obtained with:
		#   v=yade.qt.views()[0];
		#   v.lookAt, v.viewDir, v.eyePosition, v.upVector
		# see scripts/checks-and-tests/gui/testGuiMixedAlphaBox.py for example.
		self.cameraParams = cameraParams
		self.viewWaitTimeSeconds = 30.0  # sometimes the build servers are oveloaded. Let's try 30 seconds and see if it works. When testing locally put here 1 second.
		self.scrNum = 0
		# FIXME : this number 14 is hardcoded in scripts/checks-and-tests/gui/testGui.sh when testing if screenshots are present.
		self.maxTestNum = 14
		if (name != None):
			self.name = name
		else:
			self.name = ""

	def getTestNum(self):
		return self.maxTestNum

	def createEmptyFile(self, path):
		with open(path, 'a'):
			os.utime(path, None)

	def finish(self):
		#sys.exit(0)
		os._exit(0)
		#ip = get_ipython()
		#ip.ask_exit()
		#ip.pt_cli.exit()
		#quit()

	def makeNextScreenshot(self):
		time.sleep(0.25)
		subprocess.run(["/usr/bin/scrot", "-z", "scr" + "_" + self.name + "_" + str(self.scrNum).zfill(2) + ".png"])
		time.sleep(0.25)

	def clickOnScreen(self, x, y, mouseButton=1):
		time.sleep(0.25)
		subprocess.run(["/usr/bin/xdotool", "mousemove", str(x), str(y), "click", str(mouseButton), "mousemove", "restore"])
		time.sleep(0.25)

	def screenshotEngine(self):
		self.scrNum += 1
		intro = "Test '" + self.name + "', stage:" + str(self.scrNum) + " iter:" + str(O.iter)
		if (self.scrNum == 1):
			self.makeNextScreenshot()
			print(intro + " moving yade.qt.Controller()")
			yade.qt.Controller()
			yade.qt.controller.setWindowTitle('GUI test: ' + self.name)
			yade.qt.controller.setGeometry(550, 20, 500, 1100)
			yade.qt.Renderer().blinkHighlight = 'NEVER'
		if (self.scrNum == 2):
			self.makeNextScreenshot()
			print(intro + " opening yade.qt.View()")
			self.clickOnScreen(634, 1054)
			# Let's wait and see if it worked
			t = 0.0
			while (len(yade.qt.views()) == 0):
				time.sleep(0.5)
				t += 0.5
				if (t > self.viewWaitTimeSeconds):
					self.createEmptyFile("screenshots/yade_qt_View_FAILED_" + self.name)
					print("*** ERROR: failed to open QT View window by clicking ***")
					self.finish()
			# OK, now we have a view to work with.
			vv = yade.qt.views()[0]
			vv.axes = True
			if (self.cameraParams == None):
				vv.lookAt = (7.978, -4.635, 8.221)
				vv.viewDir = (-0.647, 0.441, -0.620)
				vv.eyePosition = (8.626, -5.076, 8.842)
				vv.upVector = (-0.691, 0.000, 0.721)
			else:
				vv.lookAt, vv.viewDir, vv.eyePosition, vv.upVector = self.cameraParams
			if (self.recenter):
				yade.qt.center()
		if (self.scrNum == 3):
			self.makeNextScreenshot()
			print(intro + " opening yade.qt.Inspector() , setting wire=True, setting intrGeom=True")
			self.clickOnScreen(982, 60)
			success = False
			for attempt in range(0, 30):
				try:
					yade.qt.controller.inspector.setGeometry(1050, 20, 500, 1100)
					success = True
				except:
					time.sleep(1)
			if (not success):
				self.createEmptyFile("screenshots/mouse_click_the_Inspector_open_FAILED_" + self.name)
				print("\n\n*** ERROR: self.clickOnScreen failed to open Inspector window ***\n\n")
				self.makeNextScreenshot()  # a debug screenshot.
				sys.stdout.flush()
				self.finish()
			yade.qt.controller.inspector.show()
			qt.Renderer().wire = True
			qt.Renderer().intrGeom = True
		if (self.scrNum == 4):
			self.makeNextScreenshot()
			print(intro + " changing tab to bodies, setting intrPhys=True, setting bound=True")
			self.clickOnScreen(1148, 26)
			# Previously these were used. They work, but do not always trigger crashes.
			# yade.qt.controller.inspector.close()
			# yade.qt.controller.inspector.setGeometry(1050,20,500,1100)
			# yade.qt.controller.inspector.tabWidget.setCurrentIndex(2)
			# yade.qt.controller.inspector.show()
			qt.Renderer().intrPhys = True
			qt.Renderer().bound = True
		if (self.scrNum == 5):
			self.makeNextScreenshot()
			print(intro + " clicking on interaction, setting wire=False, setting intrWire=True")
			self.clickOnScreen(1494, 55)
			qt.Renderer().wire = False
			qt.Renderer().intrWire = True
		if (self.scrNum == 6):
			self.makeNextScreenshot()
			print(intro + " changing tab to interactions, setting wire=False, setting intrWire=True")
			self.clickOnScreen(1234, 26)
			qt.Renderer().wire = False
			qt.Renderer().intrWire = True
		if (self.scrNum == 7):
			self.makeNextScreenshot()
			print(intro + " changing tab to cell, setting intrPhys=False")
			self.clickOnScreen(1312, 26)
			qt.Renderer().intrPhys = False
		if (self.scrNum == 8):
			self.makeNextScreenshot()
			print(intro + " changing tab to bodies, setting intrWire=False, setting intrGeom=False, setting intrPhys=True")
			self.clickOnScreen(1148, 26)
			qt.Renderer().intrWire = False
			qt.Renderer().intrGeom = False
			qt.Renderer().intrPhys = True
		if (self.scrNum == 9):
			self.makeNextScreenshot()
			print(intro + " changing tab to display, setting intrAllWire=True")
			self.clickOnScreen(668, 26)
			# yade.qt.controller.setTabActive('display')
			qt.Renderer().intrAllWire = True
		if (self.scrNum == 10):
			self.makeNextScreenshot()
			print(intro + " changing tab to generator, setting intrGeom=True")
			self.clickOnScreen(744, 26)
			# yade.qt.controller.setTabActive('generator')
			qt.Renderer().intrGeom = True
		if (self.scrNum == 11):
			self.makeNextScreenshot()
			print(intro + " changing tab to python, setting intrAllWire=False")
			self.clickOnScreen(821, 26)
			# yade.qt.controller.setTabActive('python')
			qt.Renderer().intrAllWire = False
		if (self.scrNum == 12):
			self.makeNextScreenshot()
			print(intro + " changing tab to simulation")
			self.clickOnScreen(580, 26)
			# yade.qt.controller.setTabActive('simulation')
		if (self.scrNum == 13):
			print(intro + " (testing of matplotlib is skipped for now...)")
			self.makeNextScreenshot()
			# FIXME: I couldn't get matplotlib to draw the plot, while screenshotting is going on.
			# self.makeNextScreenshot();
			# O.pause()
			# print(intro+" opening yade.plot.plot()")
			# plot.liveInterval=0
			# plot.plot(subPlots=False)
			# fig=yade.plot.plot();
			# time.sleep(5)
			# matplotlib.pyplot.draw()
			# time.sleep(5)
			# self.makeNextScreenshot();
			# matplotlib.pyplot.close(fig)
		if (self.scrNum == 14):
			self.makeNextScreenshot()
			self.createEmptyFile("screenshots/testGui_" + self.name + "_OK_or_Skipped.txt")
			print(intro + " exiting\n\n")
			O.pause()
			vv = yade.qt.views()[0]
			vv.close()
			yade.qt.controller.inspector.close()
			yade.qt.controller.close()
			self.finish()
