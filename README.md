# RBiQuaternions
Julia package of reduced biquaternions. 

[![Build Status](https://github.com/ivanslapnicar/RBiQuaternions.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ivanslapnicar/RBiQuaternions.jl/actions/workflows/CI.yml?query=branch%3Amain)

A reduced biquaternion $a \in\mathbb{Q}_\mathbb{R}$ can be uniquely expressed as

$$
a = a_0 + a_1 \mathbf{i} + a_2 \mathbf{j} + a_3 \mathbf{k},$$

where $a_i \in \mathbb{R}$ for
$i = 0, 1, 2, 3$, and 
$$
\mathbf{i}^2 = \mathbf{k}^2 = −1,\quad \mathbf{j}^2 = 1,\quad 
\mathbf{i}\mathbf{j} = \mathbf{j}\mathbf{i} = \mathbf{k},\quad  
\mathbf{j}\mathbf{k} = \mathbf{k}\mathbf{j} = \mathbf{i},\quad 
\mathbf{k}\mathbf{i} = \mathbf{i}\mathbf{k} = −\mathbf{j}.$$

The absolute value of $a$ is $|a| = \sqrt{a^2_0 + a^2_1+a^2_2+a^2_3}$

The conjugate of $a$ is $\bar{a}=a_0 - a_1 \mathbf{i} + a_2 \mathbf{j} - a_3 \mathbf{k}$.

Two basic singular (non-invertible) reduced biquaternions are
$$
e_1=\frac{1}{2}+\frac{1}{2}\mathbf{j}\quad \textrm{and} \quad e_2=\frac{1}{2}+\frac{1}{2}\mathbf{j}.$$
