---
layout: post
title: 자동화된 데이터 과학을 위한 트리기반 파이프라인 최적화 툴의 평가
description:
date: 2020-01-13
tags: [paper,automl,study]
---

원문 : [Evaluation of a Tree-based Pipeline Optimization Tool for Automating Data Science](https://arxiv.org/abs/1603.06212) (2016년)

### 나의 생각
이 논문은 파이썬 패키지 [TPOT (Tree-based Pipeline Optimization Tool)](https://github.com/EpistasisLab/tpot) 의 기반 논문입니다. TPOT 은 자동화된 모델 선택과 모수 최적화 등을 해주는 AutoML 파이프라인 툴이라고 할 수 있습니다. 유사한 패키지로 [auto-sklearn](https://automl.github.io/auto-sklearn/master/) 이 있으나, 몇 가지 차이점이 있습니다:
  * TPOT은 진화 알고리즘을, *auto-sklearn* 은 베이지안을 사용합니다.
  * [TPOT은 회귀에, *auto-sklearn* 은 분류에 강하다](https://medium.com/georgian-impact-blog/choosing-the-best-automl-framework-4f2a90cb1826) 는 보고가 있습니다.

저는 개인적으로 TPOT 쪽의 손을 들어주고 싶은데, 가장 큰 이유는 *[Dask](https://dask.org) 클러스터를 이용한 대량의 처리가 가능* 하기 때문입니다. 앞서, 자동화된 피처 합성을 해주는 [Featuretools](https://www.featuretools.com) 의 논문인 [심층 피처 합성 (Deep Feature Synthesis)](https://haje01.github.io/2019/12/27/deep-feature-synthesis.html) 을 소개드렸는데요, Featuretools 로 생성된 피처를 이용해 TPOT 에서 최적의 파이프라인을 구축하는 식으로 사용하면 좋을 듯 합니다.

### 주의
- 이 글은 논문의 핵심만 간추린 요약본입니다.
- 누락이나 오류가 있을 수 있으니, 정확한 내용은 꼭 원문을 참고하시기 바랍니다.

---

## 도입

* 최근 데이터와 분석 수요의 급격한 증가에 따라, 비전문가의 분석을 지원하기 위해 확장 가능하고 융통성있는 툴이 필요
* 기존 ML 툴은 최적의 데이터 파이프라인을 찾기 위해, 많은 사전 지식과 시행을 필요로 하기에 비용이 큼

전형적인 데이터 과학자는 그림 1 처럼 머신러닝 문제에 접근

![그림 1](/assets/2020-01-13-12-15-19.png)

> 그림 1 : 일반적인 지도학습 머신러닝 프로세스의 묘사.
>
> 데이터의 모델을 찾기 전 :
> - 초기 탐색적 분석 : 누락되거나 레이블이 잘못 지정된 데이터 찾기 등을 통해 데이터를 준비
> - 데이터 정리 (Data Cleaning) : 문제가 있는 기록을 수정하거나 제거
> - 피처 전처리 (Feature Preprocessing) : 피처 정규화 등
> - 피처 선택 (Feature Selection) : 모델링에 유용하지 않은 피처를 제거
> - 피처 구축 (Feature Construction) : 기존 데이터에서 새 피처를 만듦
>
> 등의 방법으로 데이터를 모델링에 보다 적합한 형식으로 변형한다. 그다음 :
> - 모델 선택 (Model Selection) - 데이터에 맞는 머신러닝 모델을 선택
> - 모수 최적화 (Parameter Optimization) : 모델이 데이터에서 가장 정확한 분류를 수행할 수 있도록 모델 모수를 선택
> - 모델 유효성 검증 (Model Validation) : 제외된 홀드 아웃 데이터 셋에서 모델의 성능을 테스트하여, 모델 예측이 수정되지 않은 데이터 셋에 일반화되도록 검증
>
> 을 수행한다. 파이프 라인에서 회색 영역이 TPOT 에 의해 자동화되는 파이프 라인의 단계를 나타낸다.

최근 *진화 알고리즘 (Evolutionary Algorithms)* 은 다양한 문제 풀이에서 인간을 능가해 왔다.

* 진화 알고리즘을 통한 머신러닝 파이프라인 설계 자동화가 가능할까?
* 우리 연구는 진화 알고리즘에 기반한 *트리 기반 파이프라인 최적화 툴 (Tree-based Pipeline Optimization Tool, TPOT)* 에 관한 것
* TPOT 은 진화 알고리즘 중 *유전자 프로그래밍 (Genetic Programming)* 을 이용해, 데이터 변환 및 머신러닝 모델을 최적화


## 관련 연구
- 지금까지의 머신러닝 자동화는 파이프라인 일부의 최적화에 집중
  - 초모수 최적화 (Hyper-parameter Optimization)
  - 피처 구축 (Feature Construction) : [심층 피처 합성](https://haje01.github.io/2019/12/27/deep-feature-synthesis.html) 을 활용한 *데이터 과학 머신*

- 최근, Fuerer 등은 베이지안 최적화를 사용하여 피처 전처리, 모델 및 모델 모수의 최적 조합을 찾아주는 *auto-sklearn* 을 개발
  - 그러나, 하나의 데이터 전처리, 피처 전처리 그리고 모델만을 포함하는 파이프라인 셋을 탐색
  - 따라서, 임의의 큰 파이프라인을 찾을 수 없음


## 구현 방법

### 파이프라인 연산자

TPOT 에 구현된 네 가지 *파이프라인 연산자 (Pipeline Operator)* 를 소개. 모든 연산자는 [scikit-learn](https://scikit-learn.org/stable/) 을 사용해 구현

* 전처리기 (Preprocessors) : `StandardScaler`, `RobustScaler`, `PolynomialFeatures`

* 분해 (Decomposition) : `RandomizedPCA`

* 피처 선택 (Feature Selection) : `Recursive Feature Elimination (RFE)`, `SelectKBest`, `SelectPercentile`, `VarianceThreshold`

* 모델 (Model) : `DecisionTree`, `RandomForestClassifier`, `GradientBoostingClassifier`, `SVM`, `LogisticRegression`, `KNeighbor`

### 트리 기반 파이프라인의 조립

* 모든 연산자들을 결합하는 유연한 파이프라인 구조를 위해 그림 2 처럼 파이프라인을 *트리 (Tree)* 로 구현

![그림 2](/assets/2020-01-13-13-58-12.png)

> 그림 2. 트리 기반 머신러닝 파이프라인의 예. 데이터 셋은 피처를 연속적으로 더하고, 제거하고, 수정하는 파이프라인 연산자를 통과한다. 결합 연산자 (Combine Operator) 는 개별 데이터 카피본을 결합하고, 최종 분류를 위해 분류기에 제공

* 트리 기반 파이프라인은 하나 이상 입력 데이터 셋 복사본을 트리의 리프로 가짐
* 그것을 네 가지 연산자 (전처리, 분해, 피처 선택, 모델링) 중 하나에 건넴
* 하나 이상의 데이터 셋의 복사본이 있으면, 결합 연산자를 통해 하나로 합침

* 데이터 셋이 모델링 연산자를 통해 통과될 때, 이전 예측을 덮어 쓰도록 최신의 분류 결과가 저장
* 데이터 셋이 파이프라인을 통해 완전히 처리되면, 최종 예측이 예측 성능 평가에 사용
* 데이터는 *층화된(Stratified)* 훈련 셋 75% 와, 테스트 셋 25% 로 나눔

* 트리기반 파이프라인은 *임의의 파이프라인 표현* 이 가능

### 진화하는 트리기반 파이프라인

* 트리기반 파이프라인 생성 및 최적화에, 파이썬 패키지  [DEAP](https://deap.readthedocs.io/en/master/) 으로 구현된 *유전 프로그래밍 (Genetic Programming)* 을 사용
* 파이프라인 연산자 시퀀스와 각 연산자의 모수를 모두 GP 로 진화
* 표 1 의 설정으로 표준 GP 과정을 따름

![표 1](/assets/2020-01-13-17-31-14.png)

* TPOT 파이프라인의 평가는 분류 정확도 (Accuracy) 로
* *TPOT-Pareto* 는 [파레토 최적화](https://ko.wikipedia.org/wiki/파레토_최적) 를 사용하는 변종
  * 최종 정확도와 파이프라인의 복잡도를 동시에 최적화 (예 : 연산자의 수)
* GP 는 진화의 세대마다 파이프라인을 조금씩 수정
  * 최종적으로 가장 우수한 파이프라인을 선택

### GAMETES 시뮬레이션 데이터 셋

* GAMETES 는 순수하고 엄격한 *상위 (Epistatic)* 유전자 모델을 생성하기 위한 오픈소스 패키지
* 이를 이용해 12 가지 유전자 모델 및 360 개 연관 데이터 셋으로 TPOT 을 평가

### UCI 벤치마크 데이터 셋

* 잘 알려진 [UC-Irvine Machine Learning Repository](https://archive.ics.uci.edu/ml/index.php) 의 9 가지 지도학습 데이터 셋으로 TPOT 을 평가
* 다양한 영역에서 TPOT 의 성능을 보여주기 위한 것
    * *Hill-Valley* : 지형 데이터에서 언덕과 계곡을 예측
    * *breast-cancder-wisconsin* : 유방암 여부를 예측
    * *car-evaluation* : 네 가지 자동차 구매 허용도를 예측
    * *glass* : 다양한 유리 데이터에서 종류를 예측
    * *ionosphere* : 고주파 안테나의 불량을 예측
    * *spambase* : 이메일 단어 빈도에서 스팸 여부 예측
    * *wine-quality* : 적/백 포도주의 11 단계 품질 예측

## 결과

TPOT 의 분류 성능을 다양한 컨트롤과 비교

* 첫 번째는 500 개 결정 트리를 가진 랜덤 포리스트
* 두 번째는 *TPOT-Random*
  * 같은 수의 파이프라인이 *임의* 로 생성됨
  * 유도된 (guided, GP로 얻은) 검색이 파이프라인 최적화에 유용한지 확인하기 위함
* 추가적으로 *TPOT-Pareto* 와 비교
  * 우수하면서 가장 작은 파이프라인을 찾기위해 *파레토 최적화* 를 이용

그림 3 은 GAMETES 데이터 셋에서 네 가지 실험을 비교

![그림 3](/assets/2020-01-14-10-39-51.png)

> 그림 3. 다양한 데이터 셋 크기와 난이도에서 TPOT 의 성능 비교
>
> * 각 서브 플롯은 테스트 데이터에서 교차 검증 정확도의 분포를 표시
> * 각 노치 박스는 30 가지 데이터 셋의 표본에 해당 (노치는 중간값의 95% 신뢰 구간)
> * 4 종류 비교
>   * *Random Forest* : 500 개 결정 트리의 랜덤 포리스트
>   * *TPOT (random search)* : 임의 파이프라인 생성
>   * *TPOT* : 유도된 검색 TPOT
>   * *TPOT (Pareto)* : 파레토 최적화된 TPOT
>
> * 각 서브 플롯은 다양한 GAMETES 설정을 표시
>   * x 축은 데이터 셋의 레코드 수
>   * y 축은 모델의 유전 가능성 (heritability, 높은 유전 가능성은 노이즈를 줄임)
> * 우상단 (높은 유전 가능성 모델에서 많은 데이터 셋이 생성) 은 쉽고 , 좌하단 (낮은 유전 가능성 모델에서 작은 데이터 셋이 생성) 으로 갈수록 어려움

* RF 는 가장 쉬운 경우 (샘플 크기: 1660, 유전 가능성: 0.4) 에서도 특성간 상위 상호 작용 발견에 실패했고, 평균 63% 정확도
* 반면, 모든 버전의 TPOT 은 80% 이상의 정확도
  * 피처 전처리와 모델링을 통해 상위 상호 작용 발견에 성공한 것
* *TPOT-Pareto* 는 정확도 분포의 분산이 낮음
  * 일관성있게 효율적인 파이프라인을 발견

그림 4 는 UCI 벤치마크 데이터 셋에 대한 네 가지 실험을 비교

![그림 4](/assets/2020-01-14-11-34-20.png)

> 그림 4. UCI 데이터 셋에 대한 TPOT 성능 비교.
>
> * 각 서브 플롯은 테스트 셋에 대한 교차 검증 정확도의 분포를 표시
> * 각 노치 박스는 데이터 셋에서 30 가지 교차 검증 분야의 샘플을 표시
> * *TPOT (random search)* 플롯의 일부가 없는 것은 120 시간 내 성공하지 못했기 때문

* TPOT 은 대부분 데이터 셋에서 RF 수준의 정확도를 달성.
  * 특히 *Hill-Valley* 와 *car-evaluation* 에서 높은 정확도를 얻음.
* *TPOT-Pareto* 는 모든 30 가지 분야에서 100% 정확도를 달성해, TPOT 을 능가
* GAMETES 와 비슷하게 *TPOT-Random* 은 유도된 TPOT 과 비슷한 성능
* 그러나, 임의 파이프라인 생성은 불필요하게 복잡한 파이프라인 때문에 유도된 검색보다 훨씬 느림
  * *Hill-Valley* 와 *spambase* 를 120 시간 내에 마치지 못함
* *TPOT-Pareto* 가 가장 컴팩트하면서 괜찮은 정확도의 파이프라인을 생성함

## 논의

* 자동화된 머신러닝 파이프라인 설계 및 최적화가 사전 지식 없이도 기초적 머신 러닝에 비해 큰 향상을 제공할 수 있다는 것을 보임
* 그러나, 자동화된 파이프라인이 데이터 과학자나 머신러닝 실무자를 완전히 대체할 수는 없음
* 우리는 TPOT 이 데이터를 탐색하고, 새 피처를 찾고, 파이프라인을 제안하는 *"데이터 과학 조수"* 가 되기를 바람
* *TPOT-Pareto* 가 생성하는 최적화된 파이프라인은, 이해하기 용이하고 실제 서비스에 적합
* 우리는 가까운 미래에 *auto-sklearn* 또는 다른 휴리스틱을 이용해 유망한 파이프라인의 모집단을 구성하도록 시도할 계획

---
