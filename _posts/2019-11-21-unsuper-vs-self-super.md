---
layout: post
title: 자율(Unsupervised)과 자기지도(Self-Supervised) 학습의 차이?
description:
date: 2019-11-21
tags: [dl, faq]
---

자율 학습은 전체 데이터셋을 학습에 사용하고, 자기지도 학습은 데이터 일부에서 나머지를 예측하도록 학습한다.

* 지도 학습은 전체 데이터 셋의 원래 형태에서 **구조**(클러스터, 밀도, 잠재 표현 등)를 학습하려 한다.
* 자기지도 학습은 데이터 저차원의 **역학**을 학습하려 한다. (예: 그레이 이미지에서 색을 예측)

<https://www.quora.com/What-is-the-difference-between-self-supervised-and-unsupervised-learning>