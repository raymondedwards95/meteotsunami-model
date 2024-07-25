# Simulations for thesis

This part is aimed at simulating meteotsunamis along the Dutch coast.
The goal is to investigate how different forcings affect the generated waves.
Here the physical parameters are important, for example bottom topography and speed and size of the low pressure system.

Two different types of pressure disturbances are used. The first type is centered around a point depending on `y`, `t` and `x0`, while the second type is centered around a line extending in the `x`-direction depending on `y` and `t`.

### Parameters 1

Pressure distribution 1 is given by
$ p \left( x, y, t \right) = p_0 \; \left( 1 - \exp\left( - \frac{t}{t_0} \right) \right) \; \exp\left( - \frac{(x - x_0)^2 + (y - U t)^2}{a^2} \right) $.

| Number | Speed `U` | Radius `a` | Shift `x0` |
| :----- | :-------- | :--------- | :--------- |
| 00     | 5         | 10000      | 0          |
| 01     | 15        | 10000      | 0          |
| 02     | 25        | 10000      | 0          |
| 03     | 35        | 10000      | 0          |
| 04     | 45        | 10000      | 0          |
| 05     | 55        | 10000      | 0          |
| 06     | 5         | 20000      | 0          |
| 07     | 15        | 20000      | 0          |
| 08     | 25        | 20000      | 0          |
| 09     | 35        | 20000      | 0          |
| 10     | 45        | 20000      | 0          |
| 11     | 55        | 20000      | 0          |
| 12     | 5         | 30000      | 0          |
| 13     | 15        | 30000      | 0          |
| 14     | 25        | 30000      | 0          |
| 15     | 35        | 30000      | 0          |
| 16     | 45        | 30000      | 0          |
| 17     | 55        | 30000      | 0          |
| 18     | 5         | 10000      | 50000      |
| 19     | 15        | 10000      | 50000      |
| 20     | 25        | 10000      | 50000      |
| 21     | 35        | 10000      | 50000      |
| 22     | 45        | 10000      | 50000      |
| 23     | 55        | 10000      | 50000      |

### Parameters 2

Pressure distribution 2 is given by
$ p \left( x, y, t \right) = p_0 \; \left( 1 - \exp\left( - \frac{t}{t_0} \right) \right) \; \exp\left( - \frac{(y - U t)^2}{a^2} \right) $.

| Number | Speed `U` | Radius `a` | Shift `x0` |
| :----- | :-------- | :--------- | :--------- |
| 50     | 5         | 10000      | 0          |
| 51     | 15        | 10000      | 0          |
| 52     | 25        | 10000      | 0          |
| 53     | 35        | 10000      | 0          |
| 54     | 45        | 10000      | 0          |
| 55     | 55        | 10000      | 0          |
| 56     | 5         | 20000      | 0          |
| 57     | 15        | 20000      | 0          |
| 58     | 25        | 20000      | 0          |
| 59     | 35        | 20000      | 0          |
| 60     | 45        | 20000      | 0          |
| 61     | 55        | 20000      | 0          |
| 62     | 5         | 30000      | 0          |
| 63     | 15        | 30000      | 0          |
| 64     | 25        | 30000      | 0          |
| 65     | 35        | 30000      | 0          |
| 66     | 45        | 30000      | 0          |
| 67     | 55        | 30000      | 0          |
