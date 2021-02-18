#-*- coding: utf-8 -*-

# 1. 用于下载原数据
import urllib

# 2. 使用shapes可以画出复杂的形状
from urllib.request import urlopen

from reportlab.graphics.renderPDF import draw
from reportlab.graphics.shapes import *

# 3. chars包里包含许多常用的图形
from reportlab.graphics.charts.lineplots import LinePlot, ScatterPlot

# 4. 用于渲染PDF文件
from reportlab.graphics import renderPDF
from reportlab.pdfgen.canvas import Canvas

c = Canvas("report3.pdf")
# 5. 初始化坐标原点
for i in [1,2,3]:
    drawing = Drawing(400, 200)

    # 7. 提取用于画图的数据
    pred = [i for i in range(0,10)]
    high = [i*i for i in range(0,10)]
    low = [i+5 for i in range(0,10)]
    times = [i for i in range(0,10)]

    # 8. zip()是Python的一个内建函数，它接受一系列可迭代的对象作为参数，将对象中对应的元素打包成一个个tuple（元组），然后返回由这些tuples组成的list（列表）。若传入参数的长度不等，则返回list的长度和参数中长度最短的对象相同。

    lp = LinePlot()
    lp.x = 10
    lp.y = 10
    lp.height = 400
    lp.width = 200
    # lp.lineLabelFormat = '%2.0f'
    # lp.xValueAxis.labelTextFormat = '%2.1f'
    # lp.xLabel='id'
    # lp.yLabel='num'
    lp.xValueAxis._valueMin=0
    lp.xValueAxis._valueMax=20
    lp.yValueAxis._valueMin = 0
    lp.yValueAxis._valueMax=100
    lp.xValueAxis.valueSteps=[i for i in range(0,10)]
    lp.data = [list(zip(times, pred)), list(zip(times, high)), list(zip(times, low))]
    lp.lines[0].strokeColor = colors.blue
    lp.lines[1].strokeColor = colors.red

    lp.lines[2].strokeColor = colors.green
    drawing.add(String(250, 150, 'x:id,y:phenoType', fontSize=14, fillColor=colors.red))

    drawing.add(lp)

    #renderPDF.drawToFile(drawing, 'report3.pdf')


    # drawing.add(String(250, 150, 'Sunspots', fontSize=14, fillColor=colors.red))
    renderPDF.draw(drawing, c, 100, 100, showBoundary=False)
    c.showPage()
#
#     #draw(drawing, c,400,200)
c.save()