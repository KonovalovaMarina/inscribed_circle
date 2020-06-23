import typing as t
import math
import sys
from PyQt5.QtWidgets import QWidget, QApplication
from PyQt5.QtGui import QPainter

EPS = pow(10, -9)
INF = sys.maxsize
steps = 60

# нахождение радиуса и центра вписанной окружности


class Point:
    def __init__(self, x: float, y: float) -> None:
        self.x = x
        self.y = y


class Line:
    def __init__(self, a, b, c) -> None:
        self.a: float = a
        self.b: float = b
        self.c: float = c


def distance(x: float, y: float, l: Line) -> float:
    return abs(x * l.a + y * l.b + l.c)


def radius(x: float, y: float, l: t.List[Line]) -> float:
    n: int = len(l)
    res: float = INF
    for i in range(n):
        res = min(res, distance(x, y, l[i]))
    return res


def y_radius(x: float, a: t.List[Point], l: t.List[Line]) -> t.Tuple[float, float]:
    n: int = len(a)
    ly: float = INF
    ry: float = -INF
    for i in range(n):
        x1: float = a[i].x
        x2: float = a[(i + 1) % n].x
        y1: float = a[i].y
        y2: float = a[(i + 1) % n].y
        if x1 == x2:
            continue
        if x1 > x2:
            x1, x2 = x2, x1
            y1, y2 = y2, y1
        if x1 <= x + EPS and x - EPS <= x2:
            y: float = y1 + (x - x1) * (y2 - y1) / (x2 - x1)
            ly = min(ly, y)
            ry = max(ry, y)
    for sy in range(steps):
        diff: float = (ry - ly) / 3
        y1: float = ly + diff
        y2: float = ry - diff
        f1: float = radius(x, y1, l)
        f2: float = radius(x, y2, l)
        if f1 < f2:
            ly = y1
        else:
            ry = y2

    return radius(x, ly, l), ly


def main():
    n: int = int(input('Введите количество вершин, а затем сами вершины в порядке обхода по многоугольнику: '))
    #n = 5
    #a = [Point(-2, -2), Point(2, -2), Point(3, 1), Point(0, 3), Point(-3, 1)]
    temp = [input().split() for _ in range(n)]
    a: t.List[Point] = [Point(int(x), int(y)) for x, y in temp]
    l: t.List[Line] = [Line(0, 0, 0) for _ in range(n)]
    for i in range(n):
        l[i].a = a[i].y - a[(i + 1) % n].y
        l[i].b = a[(i + 1) % n].x - a[i].x
        sq: float = math.sqrt(l[i].a * l[i].a + l[i].b * l[i].b)
        l[i].a /= sq
        l[i].b /= sq
        l[i].c = - (l[i].a * a[i].x + l[i].b * a[i].y)

    lx: float = INF
    rx: float = -INF
    for i in range(n):
        lx = min(lx, a[i].x)
        rx = max(rx, a[i].x)

    for sx in range(steps):
        diff: float = (rx - lx) / 3
        x1: float = lx + diff
        x2: float = rx - diff
        f1, _ = y_radius(x1, a, l)
        f2, _ = y_radius(x2, a, l)
        if f1 < f2:
            lx = x1
        else:
            rx = x2
    ans, ly = y_radius(lx, a, l)
    return ans, lx, ly, a


class Example(QWidget):

    def __init__(self):
        super().__init__()
        self.initUI()
        self.frame = 10  # отступ от верхнего и нижнего краёв окна
        self.setFixedSize(800, 800)
        self.ymin = -10
        self.ymax = 10
        self.xmin = -10
        self.xmax = 10
        self.maxy = self.size().height() - self.frame
        self.maxx = self.size().width() - self.frame

    def initUI(self):
        self.setWindowTitle('Function')
        self.show()

    def paintEvent(self, e):
        qp = QPainter()
        qp.begin(self)
        self.draw_function(qp)
        qp.end()

    # рисование окружности алгоритмом Брезенхема

    def draw_circle(self, qp, x, y, r):
        disp_x = x
        disp_y = y
        x = 0
        y = r
        delta = (1 - 2 * r)
        error = 0
        while y >= 0:
            qp.drawPoint(disp_x + x, disp_y + y)
            qp.drawPoint(disp_x + x, disp_y - y)
            qp.drawPoint(disp_x - x, disp_y + y)
            qp.drawPoint(disp_x - x, disp_y - y)

            error = 2 * (delta + y) - 1
            if (delta < 0) and (error <= 0):
                x += 1
                delta = delta + (2 * x + 1)
                continue
            error = 2 * (delta - x) - 1
            if (delta > 0) and (error > 0):
                y -= 1
                delta = delta + (1 - 2 * y)
                continue
            x += 1
            delta = delta + (2 * (x - y))
            y -= 1

    # преобразование координат в экранные

    def translate_x(self, point_x) -> int:
        return (-point_x - self.xmax) * self.maxx / (self.xmin - self.xmax) + self.frame / 2

    def translate_y(self, point_y) -> int:
        return (point_y - self.ymax) * self.maxy / (self.ymin - self.ymax) + self.frame / 2

    def draw_coordinate_axis(self, qp):
        # оси
        qp.drawLine(0, self.translate_y(0), self.maxx, self.translate_y(0))
        qp.drawLine(-self.xmin * self.maxx / (self.xmax - self.xmin) + self.frame / 2, 0,
                    -self.xmin * self.maxx / (self.xmax - self.xmin) + self.frame / 2, self.maxy + self.frame)

        # единичные отрезки по Х
        for i in range(self.xmin, self.xmax):
            qp.drawLine((i - self.xmin) * self.maxx / (self.xmax - self.xmin) + self.frame / 2,
                        (0 - self.ymax) * self.maxy / (self.ymin - self.ymax) + self.frame / 2 - 2,
                        (i - self.xmin) * self.maxx / (self.xmax - self.xmin) + self.frame / 2,
                        (0 - self.ymax) * self.maxy / (self.ymin - self.ymax) + self.frame / 2 + 2)

        # единичные отрезки по Y
        for i in range(int(self.ymin) - 1, int(self.ymax) + 2):
            qp.drawLine(-self.xmin * self.maxx / (self.xmax - self.xmin) - 2 + self.frame / 2,
                        (i - self.ymax) * self.maxy / (self.ymin - self.ymax) + self.frame / 2,
                        -self.xmin * self.maxx / (self.xmax - self.xmin) + 2 + self.frame / 2,
                        (i - self.ymax) * self.maxy / (self.ymin - self.ymax) + self.frame / 2)

    def draw_function(self, qp):
        radius, center_x, center_y, points = main()

        # оси координат
        self.draw_coordinate_axis(qp)

        # выпуклый многоугольник
        for i in range(-1, len(points)-1):
            qp.drawLine(self.translate_x(points[i].x), self.translate_y(points[i].y),
                        self.translate_x(points[i+1].x), self.translate_y(points[i+1].y))

        self.draw_circle(qp, self.translate_x(center_x), self.translate_y(center_y),
                         radius*((self.size().height()-self.frame)/(self.ymax-self.ymin)))


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = Example()
    sys.exit(app.exec_())
