# Моделирование концентрационных полей при утечке в замкнутом помещении

Этот проект моделирует распространение вещества от локальной утечки в 2D-прямоугольной области и показания течеискателя, движущегося по траектории сканирования.

## 1. Основное уравнение модели

Используется уравнение адвекции-диффузии с источником:

```math
\frac{\partial C}{\partial t} = D\left(\frac{\partial^2 C}{\partial x^2} + \frac{\partial^2 C}{\partial y^2}\right) - u_x\frac{\partial C}{\partial x} - u_y\frac{\partial C}{\partial y} + Q\,S(x,y)
```

где:

- $C(x,y,t)$ — концентрация вещества,
- $D$ — коэффициент диффузии,
- $u_x, u_y$ — компоненты скорости потока воздуха,
- $Q$ — мощность источника утечки,
- $S(x,y)$ — пространственное распределение источника.

## 2. Источник утечки

Источник задан гауссовым распределением:

```math
S(x,y) = \exp\left(-\frac{(x-x_{\text{leak}})^2 + (y-y_{\text{leak}})^2}{2\sigma^2}\right)
```

Далее выполняется нормировка:

```math
S \leftarrow \frac{S}{\sum\limits_{i,j} S_{i,j}}
```

Это делает интегральную подачу вещества менее чувствительной к размеру расчетной сетки.

## 3. Пространственная дискретизация

Для внутренних узлов применяются центральные разности:

```math
\frac{\partial C}{\partial x}\Big|_{i,j} \approx \frac{C_{i,j+1}-C_{i,j-1}}{2\,\Delta x},
\qquad
\frac{\partial C}{\partial y}\Big|_{i,j} \approx \frac{C_{i+1,j}-C_{i-1,j}}{2\,\Delta y}
```

```math
\frac{\partial^2 C}{\partial x^2}\Big|_{i,j} \approx \frac{C_{i,j+1}-2C_{i,j}+C_{i,j-1}}{\Delta x^2},
\qquad
\frac{\partial^2 C}{\partial y^2}\Big|_{i,j} \approx \frac{C_{i+1,j}-2C_{i,j}+C_{i-1,j}}{\Delta y^2}
```

## 4. Шаг по времени (явная схема Эйлера)

Обновление концентрации на шаге $n \to n+1$:

```math
C^{n+1} = C^n + \Delta t\left[D\left(\frac{\partial^2 C}{\partial x^2}+\frac{\partial^2 C}{\partial y^2}\right)-u_x\frac{\partial C}{\partial x}-u_y\frac{\partial C}{\partial y}+Q\,S\right]
```

## 5. Ограничение шага по времени (устойчивость)

Используются оценки:

```math
\Delta t_{\text{diff}} = 0.25\,\frac{\min(\Delta x^2,\Delta y^2)}{D}
```

```math
\Delta t_{\text{adv}} = 0.45\,\frac{1}{\max\left(\frac{|u_x|}{\Delta x},\frac{|u_y|}{\Delta y},\varepsilon\right)}
```

```math
\Delta t = \min(\Delta t_{\text{diff}},\Delta t_{\text{adv}})
```

(в коде $\varepsilon=10^{-12}$ для защиты от деления на ноль).

## 6. Граничные условия

Для непроницаемых стен принимается условие Неймана:

```math
\frac{\partial C}{\partial n}=0
```

В численной реализации это задается копированием ближайших внутренних значений на границы.

## 7. Модель показаний датчика

Концентрация в текущей точке прибора $(x_s,y_s)$ вычисляется билинейной интерполяцией по 4 соседним узлам:

```math
C_s = (1-t_x)(1-t_y)C_{00} + t_x(1-t_y)C_{10} + (1-t_x)t_yC_{01} + t_xt_yC_{11}
```

где $t_x, t_y \in [0,1]$ — локальные координаты внутри ячейки сетки.

---

Файл с реализацией модели: `leak_detector_simulation.m`.
