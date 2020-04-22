---
layout: post
title: 스네이크메이크 (Snakemake) 튜토리얼
description:
date: 2020-04-21
tags: [knowhow, draft]
image: /assets/2020-04-20-18-21-14.png
---


데이터 과학이나 엔지니어링을 하다보면 점점 더 복잡해지는 데이터와 코드의 숲에서 헤메는 경우가 생긴다. 데이터 처리를 체계적으로 하게 해주는 워크플로우 관리 시스템을 도입하는 것은 좋은 개선책이 될 수 있다.

![Snakemake](/assets/2020-04-20-18-21-14.png)

이 글에서는 파이썬용 워크플로우 관리 시스템들 중 하나인 Snakemake (스네이크메이크) 에 대해 실습을 통해 알아보겠다. Snakemake 버전 5.14 와 유닉스 계열 (Linux/macOS) OS 관점에서 설명하겠지만, 경로 구분자 등 몇몇 차이를 고려하면 윈도우에서도 적용에 무리가 없을 것이다.

## 워크플로우 관리가 필요한 이유

요즘은 데이터 과학과 인공지능이 각광받는 시대이지만, 사실 그것을 위한 준비 과정의 어려움은 간과되기 쉽다. 그 중 하나가 데이터 엔지니어링인데, 간단한 데이터 처리라면 그다지 복잡할 것이 없겠지만, 현업에서 높은 다양한 소스의 데이터를 처리 (ETL) 하고 고도화된 피처 엔지니어링을 수행하다 보면, 데이터 엔지니어링 단계에서 체계적인 작업 관리의 필요성을 느낄 수 있다. 대표적인 문제 몇 가지를 나열하면:

* 다양한 데이터 소스간 의존성 파악이 힘듦 - 작업에 필요한 데이터가 적절히 준비되어 있나?
* 반복적인 데이터 처리작업으로 인한 비효율성 - 꼭 필요하지 않은 전처리를 반복하는 경우
* 처리 코드와 의존 관계 혼재 - 의존하는 데이터들 간 처리 코드가 혼재되어 있으면 재활용성이 낮아지고 가독성이 떨어짐

이와 같은 문제를 해결하기 위해 데이터 처리를 위한 다양한 워크플로우 관리 시스템 (Workflow Management System, WMS) 이 활용되고 있다. 대부분의 워크플로우 매니저는 **유향 비순환 그래프 (Directed Acyclic Graph, DAG)** 방식으로 구현이 되어 있다. 아래에 간단한 DAG 구조의 예를 그림로 나타내었다.

![간단한 DAG](/assets/2020-04-16-13-32-38.png)

DAG 에서 각 노드는 작업을 의미하고, 자기 작업을 수행하기 위해 필요한 의존 작업의 노드에 연결되어 있다. 위 그림에서 최종 타겟 D를 수행하기 위해, B 와 C 가 필요하고, B 가 수행되기 위해서는 A 와 C 가 수행되어야 한다. 이렇게 타겟 결과물이 요청되면 그것을 생성하기 위해 필요한 의존 작업을 재귀적으로 거슬러 올라가며 호출하는 방식으로, 대부분 WMS 가 이 DAG 를 통해 의존관계를 해결해 나가는 식으로 구현되어 있다.

요즘은 [AirFlow](https://airflow.apache.org), [Luigi](https://github.com/spotify/luigi), [Snakemake](https://snakemake.readthedocs.io/en/stable/) 등 다양한 WMS 가 쓰이고 있는데, 이 글에서는 Snakemake (스네이크메이크) 에 대해서 알아본다. 내가 생각하는 Snakemake 의 장점은 아래와 같다.

* GNU Make 형식의 컴팩트한 구문
* 데이터 의존관계 파악이 쉽다.
* 처리 코드와 의존 관계 명세를 나누어 재활용성이 높다.
* 파이썬으로 구현되어 호환성이 좋다.

요즘처럼 통합개발환경 (IDE) 이 없었던 옛날에는, 프로그램을 빌드하기 위해서 여러가지 소스 코드 파일과 중간 산출물 파일들을 직접 관리해 주어야 했는데 이것은 꽤 복잡한 과정이었다. 이를 해결하기 위해서 사용된 것이 GNU Make인데, Snakemake는 그것의 철학을 이어받은 WMS 라고 할 수 있겠다.

> 의존성을 해결하며 최종 타겟을 만들어 가는 과정을 Make 의 관례에 따라 **빌드** 라고 하겠다.

Snakemake 는 이름이 `Snakefile` 인 워커플로우 파일에 데이터의 의존 관계와 실행 정보를 코드와 별도로 기술한다. 의존 관계와 처리 코드를 나누는 것에 대해 이야기하기 위해, 아래에 데이터 처리 관점에 특화된 시각화를 해보았다. DAG 를 데이터 의존 관계와 처리 코드, 그리고 생성된 데이터로 구분해 그렸다.

![데이터 의존 관계과 처리를 함께 기술](/assets/2020-04-16-13-41-39.png)

위 그림에서 네모는 처리 코드이고, 속이 빈 점과 화살표는 데이터와 의존관계를 나타낸다. 둘러싼 점선은 파일로 생각하면 되는데, 위의 경우 데이터 의존성과 처리가 하나의 파일에 함께 기술되었다. 일반적인 파이썬 코드로 데이터 처리를 구현하면 위와 같은 형태가 될 것이다. Snakemake 를 활용해 두 가지를 분리한 그림은 아래와 같이 되겠다.

![데이터 의존 관계와 처리를 별도로 기술](/assets/2020-04-16-13-37-16.png)

한 눈에도 좀 더 깔끔해진 것을 알 수 있다. 이렇게 하면 데이터 의존 관계를 파악하기 쉽고, 코드는 재활용성이 높아지는 효과를 기대할 수 있다.

Snakemake 를 익히는데 가장 좋은 것은 [Snakemake 공식 문서](https://snakemake.readthedocs.io) 라고 생각한다. 그런데 아쉬운 점은, Snakemake의 개발자가 Bio Informatics 를 전문으로 하고 있어서, 관련 예제가 바이오 쪽의 툴과 데이터를 이용하고 있다는 점이다. 아무래도 그 분야에 익숙하지 않은 사람들에게는 생소하기 때문에, 이 글에서는 가급적 일반적인 예로 설명하려 한다.

## 빠르게 시작하기

간단한 데이터 처리 예를 통해 Snakemake의 기본 개념을 확인해 보자. 먼저 아래와 같이 Snakemake 를 설치한다.

```
$ pip install snakemake
```

> Snakemake 는 파이썬 3.5 이상 버전을 필요로 한다. Python2 와 Python3를 함께 사용 중이라면 `pip3` 로 설치해야 할 것이다.

앞으로 설명할 예제는 각 파일에 들어있는 단어의 수를 세고, 최종적으로 파일별 단어 수 그래프를 만드는 것을 목표로 한다.

먼저 데이터 디렉토리가 다음과 같이 구성되어 있다고 하자:

```
data/
    A.txt
    B.txt
    C.txt
```

파일의 내용은 다음과 같다.

`A.txt`
```
Hello.
```

`B.txt`
```
This is your captain speaking.
```

`C.txt`
```
Welcome aboard.
```

먼저 워크플로우 파일인 `Snakefile` 이 필요하다. Snakemake 는 이 워크플로우 파일의 내용대로 동작하게 된다.

일반적으로 워크플로우는 하나 이상의 작업 노드들이 연결된 DAG 형식인데, Snakemake 에서는 이 노드를 **규칙 (Rule)** 이라 부르며, `Snakefile` 에서 `rule` 이라는 키워드로 선언한다. 규칙은 기본적으로 다음과 같은 구조를 가진다.

```
rule 규칙 이름:
    입력 파트:
    출력 파트:
    실행 파트:
```

파트별 구체적인 키워드는 다음과 같다.

```
rule RULE_NAME:
    input:
    output:
    run, script, shell, notebook 중 선택:
```

* `rule` 에는 규칙의 이름이 온다.
* `input` 에는 규칙의 입력이 되는 하나 이상의 파일 경로나, 코드 또는 함수가 온다. 경로가 하나 이상인 경우 `,` 를 사용하여 구분한다.
* `output` 에는 규칙의 출력이 되는 하나 이상의 파일 경로가 온다. 하나 이상인 경우 `,` 를 사용하여 구분한다.
* 실행 파트에는 다음과 같은 키워드가 온다
  * `shell` - 실행할 쉘 명령을 기술
  * `run` - 파이썬 스크립트를 직접 기술
  * `script` - 외부 파이썬 스크립트 파일 경로를 지정
  * `notebook` - Jupyter 노트북 파일 경로를 지정

> `Snakefile` 은 기본적으로 파이썬 파일과 같다. 즉, 파일 상단에서 파이썬 모듈을 임포트하거나 파이썬 명령어를 사용할 수 있다. 단, 위에서 보듯 **규칙 (Rule)** 부분은 다른 구문을 사용한다.

### 파일내 단어 수 세기

이제 다음과 같은 내용으로 `Snakefile` 파일을 현재 디렉토리에 만들어 보자:

```python
rule count:
    """파일내 단어 수 세기."""
    input:
        "data/{filename}.txt"
    output:
        "temp/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"
```

위 예는 `count` 라는 규칙을 정의한다. `data/` 디렉토리 아래의 `[파일명].txt` 에 매치되는 모든 파일을 입력으로 하여 `shell` 에 기술된 쉘 명령을 수행한 후, 그 결과를 `temp/wc_[파일명].txt` 파일에 출력한다는 뜻이다.

> `wc -w` 는 파일에서 단어를 세는 Unix 계열 쉘 명령어다. 예를 들면 다음과 같다.
> ```
> $ wc -w data/B.txt
>
>        5 data/B.txt
> ```
> `B.txt` 에 다섯 단어가 있음을 알 수 있다.

입력 및 출력에서 `{filename}` 부분은 **와일드카드** 로 불리는데, 패턴이 매칭되는 모든 파일에서 대체되는 변수이다.

Snakemake 의 실행은 `Snakefile` 이 있는 디렉토리에서 다음과 같이 한다.

```
$ snakemake [빌드 타겟] -j [숫자]
```

`-j` 옵션은 `--jobs` 또는 `--cores` 의 단축형으로 숫자 인자를 받는다. 그것은 Snakemake 빌드에서 최대 몇 개의 CPU 코어를 활용할 것인지 명시하는 것이다. 많은 코어를 사용하면 DAG 에 따라 의존성이 없는 규칙들을 병렬로 처리하여 속도를 올릴 수 있다. 숫자를 생략하면 가용한 최대 코어를 이용한다. 빌드 타겟은 Snakemake 를 실행하여 얻고자 하는 최종 타겟 파일을 지정한다.

예제를 위해 다음과 같이 실행해보자:

```
$ snakemake temp/wc_A.txt -j
```

위 명령은 `temp/wc_A.txt` 를 타겟으로 빌드하라는 것인데, 이를 위해서는 출력 ( `output` ) 이 `temp/wc_A.txt` 에 매칭되는 규칙이 필요하다. 여기서 그 규칙은 `count` 이다. 이 규칙에서 출력이 매칭되면 와일드카드의 `{filename}` 은 `A` 로 배정되며, 이에 의해 입력 ( `input` ) 은 `data/A.txt` 로 결정된다. 결과적으로 실행은 쉘명령어 `wc -w data/A.txt > temp/wc_A.txt` 를 수행하게 된다. 실행 후 다음처럼 결과를 확인해보자:

```
$ cat temp/wc_A.txt

       1 data/A.txt
```

`data/A.txt` 파일의 단어 수가 저장되었다.

이처럼 Snakemake는 **각 단계별 결과물을 파일로 저장하는데, 이를 통해 반복 수행을 피하고 중간 결과 확인이 용이** 해진다.

### 개별 단어 수 파일 병합

다음 단계는 각 파일별 단어 수를 하나의 `.csv` 형식 파일로 저장하는 것이다. 이를 위해 `Snakefile` 에 아래처럼 `concat` 규칙을 추가한다.

```python
rule concat:
    """개별 단어 수 파일을 병합."""
    input:
        "temp/wc_A.txt",
        "temp/wc_B.txt",
        "temp/wc_C.txt"
    output:
        "temp/wc_all.csv"
    script:
        "concat.py"
```

위 규칙은 앞에서 생성한 모든 단어 수 파일을 입력으로 한다. 즉, `count` 규칙은 입력과 출력 파일이 **일대일**의 관계였지만 `concat` 규칙은 **다대일**의 관계이다.

`input` 에 모든 단어 수 파일을 명시하고 있는데, 이렇게 해도 되지만 아래처럼 `expand` 함수를 사용하면 편리하다.

```python
rule concat:
    """개별 단어 수 파일을 병합."""
    input:
        expand("temp/wc_{filename}.txt", filename=['A', 'B', 'C'])
    output:
        "temp/wc_all.csv"
    script:
        "concat.py"
```

`expand` 는 매개변수로 받은 리스트의 각 항목으로 패턴을 확장한다. 위의 경우, 결과는 값은 아래와 같다.

```python
["temp/wc_A.txt", "temp/wc_B.txt", "temp/wc_C.txt"]
```

> `expand` 의 첫 인자 내 `{filename}` 은 와일드카드가 아니고, 매개변수 리스트의 각 항목이 대입되는 위치이다. `expand` 인자에서 와일드카드를 쓰려면 이중 중괄호 {%raw%}`{{ }}`{%endraw%} 를 사용해야 한다.

`Snakefile` 은 기본적으로 파이썬 파일이기에, 다음과 같이 상수를 정의해 이용할 수 있다.
```python
FILENAMES = ['A', 'B', 'C']

rule concat:
    """개별 단어 수 파일을 병합."""
    input:
        expand("temp/wc_{filename}.txt", filename=FILENAMES)
    output:
        "temp/wc_all.csv"
    script:
        "concat.py"
```

이 규칙에서는 세 파일을 결합하기 위해 `concat.py` 이라는 파이썬 파일을 실행해 `temp/wc_all.csv` 로 출력하는데, 그 코드는 아래와 같다.

```python
import re
import pandas as pd

# wc -w 결과 파싱용 정규식
PTRN = re.compile(r'\s*(\d+)\s[^\s]+([^\s\/]+.txt)')

# 출력용 csv 파일 오픈
with open(snakemake.output[0], 'wt') as f:
    f.write('fname, count\n')
    # snakefile에 명시된 모든 입력 파일에 대해서
    for fn in snakemake.input:
        # 파싱하고 출력
        line = open(fn, 'rt').read()
        cnt, fn = PTRN.search(line).groups()
        f.write('{}, {}\n'.format(fn, cnt))
```

Snakemake 를 통해서 실행되는 파이썬 스크립트에는 `snakemake` 라는 전역 객체가 기본적으로 생성된다. 여기에는 Snakemake 관련 실행 정보가 들어있다. 위 스크립트에서는 `Snakefile` 의 `concat` 규칙에서 정해진 입력 파일들을 얻어 오기 위해 `snakemake.input` 을, 출력 파일을 얻어 오기 위해서 `snakemake.output` 을 이용하고 있다. `Snakefile` 의 `input` 이나 `output` 에는 하나 이상의 파일이 올 수 있기에, 리스트로 취급한다.

> 위의 예처럼, 가급적 Snakemake 에서 제공하는 정보를 최대한 활용하는 것이 코드의 재활용성을 높여주는 좋은 습관이다.

이제 다음처럼 실행하면:

```
$ snakemake temp/wc_all.csv -j
```

`A.txt`, `B.txt`, `C.txt` 세 파일의 단어 수가 하나의 `temp/wc_all.csv` 파일에 저장된다. 다음처럼 확인해보자.

```
$ cat temp/wc_all.csv

fname, count
A.txt, 1
B.txt, 5
C.txt, 2
```

### 그래프 그리기

이제 위에서 생성된 `temp/wc_all.csv` 를 읽어 그래프를 그리는 마지막 규칙을 살펴보자. `Snakefile` 에 아래와 같은 규칙을 추가한다.

```python
rule plot:
    """그래프 그리기."""
    input:
        "temp/wc_all.csv"
    output:
        "temp/wc_all.png"
    script:
        "plot.py"
```

스크립트파일 `plot.py` 은 `temp/wc_all.csv` 파일을 읽어 막대 그래프를 `temp/wc_all.png` 에 저장하는데, 코드는 아래와 같다.

```python
import pandas as pd

df = pd.read_csv(snakemake.input[0], index_col='fname')
plot = df.plot(kind='bar', rot=45)
fig = plot.get_figure()
fig.tight_layout()
fig.savefig(snakemake.output[0])
```

실행하면 멋진(?) 그래프가 만들어진다.

```
$ snakemake temp/wc_all.png -j
```

![단어 그래프](/assets/2020-04-08-14-27-30.png)


### 기본 규칙 만들기

Snakemake 를 호출할 때마다 매번 타겟 파일을 지정하는 것은 번거롭다. Snakemake 는 타겟 파일이 없으면 첫 번째 규칙을 실행하는데, 다음의 규칙을 `Snakefile` 의 첫 규칙으로 추가해 기본 규칙으로 동작하게 한다.

```python
rule all:
    """기본 규칙."""
    input: "temp/wc_all.png"
```

이 규칙에는 입력만 있는데, 이 파일을 생성하는 다른 규칙을 찾는 것에서 빌드가 시작된다. 이제 아래처럼 타겟 파일없이 실행하면 자동으로 `all` 규칙이 선택되고, 그것의 입력인 `temp/wc_all.png` 를 생성하는 `plot` 규칙이 실행된다.

```
$ snakemake -j
```

지금까지 작업한 전체 `Snakefile` 의 내용은 아래와 같을 것이다.

```python
FILENAMES = ['A', 'B', 'C']

rule all:
    """기본 규칙."""
    input:
        "temp/wc_all.png"

rule count:
    """파일내 단어 수 세기."""
    input:
        "data/{filename}.txt"
    output:
        "temp/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    """개별 단어 수 파일을 병합."""
    input:
        expand('temp/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/wc_all.csv"
    script:
        "concat.py"

rule plot:
    """그래프 그리기."""
    input:
        "temp/wc_all.csv"
    output:
        "temp/wc_all.png"
    script:
        "plot.py"
```

### DAG 시각화하기

워크플로우 파일의 DAG 를 아래의 명령으로 시각화할 수 있다.

```
$ snakemake -j --dag | dot -Tpng -o dag.png
```

`dag.png` 를 보면 다음과 같을 것이다.

![DAG 시각화](/assets/2020-04-16-13-43-54.png)

### DAG 동작 돌아보기

지금까지 우리는 단계별로 진행을 했기에 중간 결과물들이 차례되로 생성되어 있었지만, 처음부터 빌드하는 경우 다음과 같이 진행된다.

1. 명시적 타겟 파일 없이 Snakemake 실행
2. 첫 규칙인 `all` 에서 `temp/wc_all.png` 파일을 입력으로 요구
3. 그것을 출력하는 `plot` 규칙에서 `temp/wc_all.csv` 파일을 입력으로 요구
4. 그것을 출력하는 `concat` 규칙에서 `temp/wc_A.txt`, `temp/wc_B.txt`, `temp/wc_C.txt` 파일을 입력으로 요구
5. 그것들을 출력하는 `count` 규칙에서 입력인 `data/A.txt`, `data/B.txt`, `data/C.txt` 의 단어를 세고 결과를 각각 `temp/wc_A.txt`, `temp/wc_B.txt`, `temp/wc_C.txt` 에 저장
6. 입력이 해결된 `concat` 규칙에서 `temp/wc_all.csv` 파일을 생성
7. 입력이 해결된 `plot` 규칙에서 `temp/wc_all.png` 파일을 생성

이러한 DAG 의존 관계를 통한 워크플로우는 순차적으로 작업을 하는 것이 아니라, 요청된 타겟 파일의 의존성을 쫓아 재귀적으로 진행되기에 일종의 **On-Demand** 또는 **Lazy Evaluation** 방식으로 볼 수 있겠다. 이 방식은 필요없는 작업이 적고, 앞 규칙의 결과 파일을 이용하기에 여러번 재시도해도 최단 시간에 결과를 확인할 수 있다.

예를 들어 앞의 `plot` 규칙에서, 좀 더 미려한 그래프를 그리기 위해서 코드를 여러번 수정해야 한다면, (Jupyter 노트북처럼 인터랙티브한 환경이 아니라면) 매번 Python 을 실행하여 결과를 확인해야 한다. 이 과정에서 앞에서 실행했던 `count` 나 `concat` 이 불필요하게 실행될 수 있다. 물론, 코드를 수정하여 각 단계에서 중간 결과를 저장하고, 중간 결과가 있는 경우 그 것을 이용하도록 할 수도 있겠으나, 아무래도 번거롭고 코드가 지저분해지게 된다. Snakemake 를 사용하면 앞 단계 규칙의 결과인 `temp/wc_all.csv`가 있으면 그래프 그리는 코드만 실행하기에 효율적이다.

## 심화

이제 Snakemake 에 대한 깊이 있는 활용 법을 알아보자.

### 두 가지 의존성

지금까지 막연하게 의존성을 이야기했으나, 사실 Snakemake 에서 다루는 의존성은 크게 두 가지로 나누어 설명할 수 있다. 그것은 **데이터 의존성** 과 **워크플로우 의존성** 이다 (이 두 용어는 내가 설명의 편의를 위해 도입한 것이다) .

#### 데이터 의존성

데이터 의존성은 타겟 파일 생성 후 입력 파일 (데이터) 이 변경되었다면 빌드가 무효화되는 것을 말한다. 변경 여부는 파일의 수정 시간을 기준으로 한다. Snakemake 는 데이터가 변하면 체크하여 자동으로 필요한 빌드를 다시 해준다. 예를 들어 아래와 같이 `temp/A.txt` 파일의 내용을 바꾸고,

```
Hello everybody.
```

Snakemake 를 다시 호출하면,

```
$ snakemake -j
```

변경된 `data/A.txt` 에 대해서만 `count` 규칙이 수행되고 `data/B.txt` 와 `data/C.txt` 에 대해서는 수행되지 않는다. 그후 갱신된 출력인 `temp/wc_A.txt` 에 의존하는 `temp/wc_all.csv` 와 그것에 의존하는 `temp/wc_all.png` 순으로 빌드가 진행되는 것을 볼 수 있다. 결과물인 그래프를 보면 `A.txt` 항목의 단어 수가 2 로 바뀐 것을 확인할 수 있다.

![단어 수 그래프 2](/assets/2020-04-09-12-57-20.png)

#### 워크플로우 의존성

워크플로우 의존성은 `Snakefile` 에 기술된 규칙이나, 그것을 처리하는 코드에 변경이 있어 빌드가 무효화되는 것을 말한다. 그런데 Snakemake 는 이런 변화는 자동으로 체크하지 못하기에, 유저가 직접 검사하여 빌드를 해주어야 한다. 다음과 같은 것들이 있다.

* 규칙의 입력 파트 변경 - `input` 의 내용이 바뀌면 무효화
  * `snakemake --list-input-changes` 로 영향 받는 출력물을 확인 가능
* 규칙의 패러미터 파트 변경 - 이후 설명할 `params` 의 내용이 바뀌면 무효화
  * `snakemake --list-params-change` 로 영향 받는 출력물을 확인 가능
* 규칙의 실행 파트 변경 - `run`, `shell`, `script` 등의 내용이 바뀌면 무효화
  * `snakemake --list-code-changes` 로 영향 받는 출력물을 확인 가능

예를 들어 다음처럼 `concat.py` 파일에 주석을 추가하고,

```python
import re
import pandas as pd

PTRN = re.compile(r'\s*(\d+)\s[^\s]+([^\s\/]+.txt)')

# 코드 갱신을 위한 주석 추가 <--
with open(snakemake.output[0], 'wt') as f:
    f.write('fname, count\n')
    for fn in snakemake.input:
        line = open(fn, 'rt').read()
        cnt, fn = PTRN.search(line).groups()
        f.write('{}, {}\n'.format(fn, cnt))

```

Snakemake 를 다시 실행해보면, 코드가 바뀌었음에도 빌드할 것이 없다고 나온다.

```
$ snakemake -j

Building DAG of jobs...
Nothing to be done.
```

이것은 Snakemake 가 코드가 바뀐 것을 인식하지 못하기 때문이다. 다음처럼 코드 변화의 영향을 명시적으로 검사할 수 있다.

```
$ snakemake --list-code-changes

Building DAG of jobs...
temp/wc_all.png
temp/wc_all.csv
```

`concat.py` 가 변했으니 그것의 출력인 `temp/wc_all.csv` 과 이 파일에 의존하는 `temp/wc_all.png` 가 모두 빌드되어야 한다고 알려준다.

코드 변화를 명시적으로 검사하고, 영향받는 타겟만 빌드를 해주기 위해서는 아래와 같이 명령하면 된다.

```
$ snakemake -j -R `snakemake --list-code-changes`
```

> `-R` 옵션은 `--forcerun` 의 단축형으로, 해당 타겟을 강제적으로 빌드하는 옵션이다.

워크플로우 내 입력 및 매개변수 변화에 대해서도 위와 같은 방식으로 대응할 수 있다.

### 사전형으로 입출력 기술하기

많은 입출력 파일을 사용하는 경우 리스트 형식으로 나열하는 것은 혼란스러울 수 있다. 이때는 `키 = 값` 형태로 입출력을 명시하면 편하다.

```
rule RULE_NAME:
    input:
        foo_in="data/foo.txt"
        boo_in="data/boo.txt"
    output:
        foo_out="temp/wc_foo.txt"
        boo_out="temp/wc_boo.txt"
```

쉘 명령에서는 `{input.foo_in}` 또는 `{output.foo_out}` 처럼, 스크립트에서는 `snakemake.input.foo_in` 또는 `snakemake.output.foo_out` 처럼 입출력을 참조한다.

### 규칙에 매개변수 이용하기

규칙에서 입출력 외에 실행 파트로 전달하고 싶은 값이 있을 때 `params` 키워드를 이용할 수 있다. `input` 처럼 하나 이상의 값이나, 코드 또는 함수가 올 수 있다. 아래는 와일드카드를 사용하는 예로, 처리된 파일 이름을 표준 출력에 표시한다.


```python
rule count:
    """파일내 단어 수 세기."""
    input:
        S3.remote("wzdat-seoul/data/{filename}.txt")
    params:
        fname="{filename}"
    output:
        "temp/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output} && echo 'Done {params.fname}'"
```

쉘 명령에서는 `params.fname` 처럼, 스크립트에서는 `snakemake.params.fname` 처럼 규칙의 매개변수를 참조한다. 실행하면 단어 수를 센 다음 다음과 같은 메시지를 출력한다.

```
Done A
```

다음처럼 `lambda` 함수를 이용할 수도 있다.

```python
rule count:
    """파일내 단어 수 세기."""
    input:
        S3.remote("wzdat-seoul/data/{filename}.txt")
    params:
        pair=lambda wildcards, output: "{} - {}".format(wildcards, output)
    output:
        "temp/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output} && echo 'Done {params.pair}'"
```

다음과 같은 메시지를 출력한다.

```
Done A - temp/wc_A.txt
```

### 스크립트 파일 합치기

앞의 예제에서는 `concat` 과 `plot` 규칙을 위해 각각 파이썬 파일을 하나씩 만들어 주었다. 만약 파이썬 파일이 너무 많아질 것 같아 하나로 통합하고 싶다면 어떻게 할까? 먼저 `Snakefile` 을 아래처럼 수정하고:

```python
FILENAMES = ['A', 'B', 'C']

rule all:
    input:
        "temp/wc_all.png"

rule count:
    input:
        "data/{filename}.txt"
    output:
        "temp/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    input:
        expand('temp/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/wc_all.csv"
    script:
        "main.py"  # <-- 스크립트 파일명 변경

rule plot:
    input:
        "temp/wc_all.csv"
    output:
        "temp/wc_all.png"
    script:
        "main.py"  # <-- 스크립트 파일명 변경

```

`concat.py` 와 `plot.py` 의 내용을 다음처럼 하나의 파이썬 파일 `main.py` 에 옮긴다.

```python
import re
import pandas as pd

PTRN = re.compile(r'\s*(\d+)\s[^\s]+([^\s\/]+.txt)')


def concat():
    with open(snakemake.output[0], 'wt') as f:
        f.write('fname, count\n')
        for fn in snakemake.input:
            line = open(fn, 'rt').read()
            cnt, fn = PTRN.search(line).groups()
            f.write('{}, {}\n'.format(fn, cnt))


def plot():
    df = pd.read_csv(snakemake.input[0], index_col='fname')
    plot = df.plot(kind='bar', rot=45)
    fig = plot.get_figure()
    fig.tight_layout()
    fig.savefig(snakemake.output[0])


if __name__ == '__main__':
    globals()[snakemake.rule]()  # <-- 현재 규칙에 맞는 함수 호출
```

파이썬 파일이 Snakemake 를 통해 호출되는 경우 `snakemake.rule` 속성에 현재 규칙의 이름이 온다는 것을 이용하였다.

### 특정 디렉토리내 모든 파일을 입력으로 하기

앞의 예에서는 어떤 입력 파일들이 있는지 알고 있다고 전제했다. `A`, `B`, `C` 가 그것이었다. 만약 확장자는 같지만 임의 파일명으로 여러 파일이 있다면 어떻게 해야할까? 파이썬의 `glob` 모듈 은 패턴에 매칭되는 파일들을 찾아주는데, Snakemake 가 제공하는 `glob_wildcards` 도 비슷한 역할을 한다.

아래의  코드를,

```python
FILENAMES = ['A', 'B', 'C']
```

다름처럼 수정해도 같다.

```python
FILENAMES, = glob_wildcards('temp/wc_{filename}.txt')
```

> 객체가 아닌 파일명 리스트를 얻기 위해 `FILENAMES` 다음에 `,` 가 있음에 주의하자.

### 와일드카드에 정규식으로 제약걸기

특정 타겟을 빌드할 때 타겟 이름을 요소별로 분리해 와일드카드에 배정하고 싶을 때가 있다. 설명을 위해, 앞의 예제를 다음과 같이 바꾸어 생각해보자:

* 일별 파일들이 `YYYYMMDD` 형식 디렉토리 아래에 저장되고,
* 일별 단어 수 그래프를 `YYYYMMDD.png` 형식으로 빌드해야 한다.

예를 들어 다음과 같은 데이터가 있다면,

```
data/
    20200401/
        A.txt
        B.txt
        C.txt

    20200402/
        A.txt
        B.txt
        C.txt

    20200403/
        A.txt
        B.txt
        C.txt
```

2020년 4월 3일의 결과물을 `snakemake temp/20200403.png -j` 명령으로 얻으려면 `Snakefile` 의 입력은 `year`, `month`, `day` 의 세 가지  와일드카드를 이용해 다음과 같이 기술될 수 있을 것이다.

```python
FILENAMES = ['A', 'B', 'C']

rule count:
    input:
        "data/{year}{month}{day}/{filename}.txt"
    output:
        "temp/{year}{month}{day}/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    input:
        expand('temp/{%raw%}{{year}}{{month}}{{day}}{%endraw%}/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/{year}{month}{day}/wc_all.csv"
    script:
        "main.py"

rule plot:
    input:
        "temp/{year}{month}{day}/wc_all.csv"
    output:
        "temp/{year}{month}{day}/wc_all.png"
    script:
        "main.py"
```

> `concat` 규칙의 `expand` 에서 와일드카드를 이용하기 위해, 앞에서 말한 이중 중괄호를 사용하였다.

그러나 Snakemake 는 `year` 가 4 자리, `month` 가 2 자리, `day` 가 2 자리인 것을 모르기 때문에, 잘못된 와일드카드 배정으로 인한 오류나 `RecursionError` 가 발생할 수도 있다. 따라서, 안전한 방법은 {%raw%}`{{와일드카드_이름, 정규식}}`{%endraw%} 형식으로 와일드카드의 패턴을 제약하는 것이다. 다음과 같이 할 수 있다.

```python
FILENAMES = ['A', 'B', 'C']

rule count:
    input:
        "data/{year}{month}{day}/{filename}.txt"
    output:
        "temp/{year,\d{4}}{month,\d{2}}{day,\d{2}}/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    input:
        expand('temp/{%raw%}{{year}}{{month}}{{day}}{%endraw%}/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/{year,\d{4}}{month,\d{2}}{day,\d{2}}/wc_all.csv"
    script:
        "main.py"

rule plot:
    input:
        "temp/{year}{month}{day}/wc_all.csv"
    output:
        "temp/{year,\d{4}}{month,\d{2}}{day,\d{2}}/wc_all.png"
    script:
        "main.py"
```

각 규칙의 출력에서 정규식으로 와일드카드의 값을 제약하고 있다. 와일드카드의 패턴을 제한하는 것은 출력 즉, `output` 에서만 사용할 수 있는데, 출력 와일드카드가 정해지면 입력은 그것을 그대로 쓰기 때문이다.

### S3에서 파일 입출력

클라우드 스토리지 서비스인 AWS S3 를 입출력 대상으로 이용할 수 있다. 앞의 예제를 S3를 이용하도록 다음처럼 바꾸어 생각해보자.

* 디렉토리에 단어 수를 셀 텍스트 파일은 `s3://my-bucket/data/` 아래에 있음
* 중간 산출물은 로컬의 `temp/` 디렉토리에 출력
* 최종 타겟은 `s3://my-bucket/result/wc_all.png` 로 저장

이를 위해 Snakemake 에서 제공하는 `S3RemoteProvider` 를 사용한다. S3 URL 에서 `s3://` 부분을 생략하고 다음처럼 수정한다.

> 예제의 `my-bucket` 은 실제 자신이 사용할 버킷 이름으로 교체해야 한다.

```python
from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
S3 = S3RemoteProvider(keep_local=True)

FILENAMES = ['A', 'B', 'C']

rule all:
    input:
        S3.remote("my-bucket/temp/wc_all.png")

rule count:
    input:
        S3.remote("my-bucket/data/{filename}.txt")
    output:
        "temp/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    input:
        expand('temp/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/wc_all.csv"
    script:
        "concat.py"

rule plot:
    input:
        "temp/wc_all.csv"
    output:
        S3.remote("my-bucket/temp/wc_all.png")
    script:
        "plot.py"
```

`S3RemoteProvider` 를 사용하면 다음과 같은 방식으로 S3 에 입출력을 할 수 있게 된다:

* 입력 파일이 S3 에 있는 경우, 먼저 Snakemake 가 동명의 로컬 디렉토리에 내려받은 뒤 스크립트는 그 파일을 이용.
* 출력 파일이 S3 에 있는 경우, 먼저 스크립트가 동명의 로컬 디렉토리에 출력한 뒤 Snakemake 는 그 파일을 S3 로 업로드 한다.

S3 입출력을 위해 임시로 사용한 로컬 파일은 더 이상 의존하는 규칙이 없으면 Snakemake 에 의해 지워진다. `S3RemoteProvier` 생성시 `keep_local=True` 으로 하면 지워지지 않는다. 빌드 수행 후 현재 디렉토리 아래 `my-bucket` 디렉토리가 만들어진 것을 발견할 수 있다. 그 내용은 다음과 같다.

```
$ ls my-bucket/
A.txt	B.txt	C.txt

$ my-bucket/temp/wc_all.png
wc_all.png
```

이렇게 로컬 디렉토리를 경유하는 방식은, 각 쉘 명령이나 스크립트가 S3 API 없이 S3 를 이용할 수 있어 편리하다. 다만 크기가 큰 파일의 경우 로컬 디스크의 용량에 주의할 필요가 있겠다.

> 위의 예는 AWS CLI 툴의 설치 및 설정이 된 것을 가정하고 있다. 만약 그렇지 않다면, `S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET")` 식으로 AWS 계정 정보를 넣어 주어야 한다.

> 다음처럼 `S3RemoteProvider` 를 통해서도 `expand` 를 수행할 수 있다.
>
> `S3.remote(expand("my-backet/data/{filename}.txt", filename=FILENAMES))`


### Jupyter 노트북 실행하기

Snakemake 는 Jupyter 노트북으로 구현된 코드를 실행 파트로 사용할 수 있다. 앞의 `concat` 규칙의 실행 파트를 `notebook` 키워드로 수정한다.

```python
"""개별 단어 수 파일을 병합."""
rule concat:
    input:
        expand('temp/wc_{filename}.txt', filename=FILENAMES)
    output:
        "temp/wc_all.csv"
    notebook:
        "concat.ipynb"
```

파이썬 파일 대신 Jupyter 노트북 파일인 `.ipynb` 가 사용되었다. 문제점은 노트북에서 개발할 때는 `snakemake` 객체가 없다는 것이다. 개발할 때는 명시적으로 파일명을 이용하고, 빌드할 때는 `snakemake` 를 참조하도록 수정할 수 있겠으나, 아무래도 실수의 여지가 많다. 하나의 아이디어는 개발시에만 동작하는 아래와 같은 막업 코드를 이용하는 것이다.

```python
# snakemake 막업 (Mockup) 만들기

class AttrDict(dict):
    def __getattr__(self, attr):
        return self[attr]
    def __setattr__(self, attr, value):
        self[attr] = value

# Snakemake 로 실행되지 않을 때만 막업 사용
if 'snakemake' not in globals():
    snakemake = AttrDict()
    # 노트북 개발시 코드 동작을 확인할 수 있는 적당한 입출력
    snakemake.input = ['temp/wc_A.txt', 'temp/wc_B.txt']
    snakemake.output = ['temp/wc_all.csv']
```

`concat` 규칙을 위한 Jupyter 노트북 파일 `concat.ipynb` 는 아래와 같은 모습이 될 것이다.

![Jupyter 노트북 코드](/assets/2020-04-22-11-25-05.png)

노트북 개발시에는 막업의 정보를, Snakemake를 통한 실행시에는 완전한 정보를 사용하게 된다.

### 다른 워크플로우 파일 사용하기

하나의 `Snakefile` 에서 다른 `Snakefile` 을 포함 (`include`) 해서 사용하는 것이 가능하다. 이 기능은 빌드를 계층화 하거나, 여러 사용자가 협업할 때 유용할 것이다. 앞의 예를 데이터 엔지니어 (u1) 와 데이터 분석가 (u2) 가 나누어 작업하는 경우로 바꾸어 생각해보자. 엔지니어는 데이터 파일을 읽어 단어수 파일 `wc_all.csv` 를 만들고, 분석가가 이를 이용해 시각화 파일 `wc_all.png` 을 만드는 식으로 가정하고, git 등으로 코드 관리에 용이하게 다음처럼 디렉토리를 구성한다.

```
data/
    A.txt
    B.txt
    C.txt
temp/
u1/
    Snakefile
    concat.py
u2/
    Snakefile
    plot.py
```

여기서 `data` 및 `temp` 디렉토리는 공유하는 것으로 한다. 먼저, 데이터 분석가 (u1) 의 `Snakefile` 은 다음과 같다.

```python
FILENAMES = ['A', 'B', 'C']

rule u1all:
    """사용자 1 의 기본 규칙."""
    input:
        "../temp/wc_all.csv"

rule count:
    """파일내 단어 수 세기."""
    input:
        lambda x: "../data/{filename}.txt"
    output:
        "../temp/wc_{filename}.txt"
    shell:
        "wc -w {input} > {output}"

rule concat:
    """개별 단어 수 파일을 병합."""
    input:
        expand('../temp/wc_{filename}.txt', filename=FILENAMES)
    output:
        "../temp/wc_all.csv"
    script:
        "concat.py"
```

앞에서 본 예와 크게 다르지 않으나, `concat` 규칙까지만 구현하고 기본 규칙의 이름을 `u1all` 로 한 것이 눈에 띄인다. 이는 다른 워크플로우 파일을 포함해 사용할 때, 같은 이름의 규칙이 있으면 충돌이 일어나기 때문이다.

데이터 분석가 (u2) 의 `Snakefile` 은 아래와 같다.

```python
include: "../u1/Snakefile"

rule u2all:
    """사용자 2 의 기본 규칙."""
    input:
        "../temp/wc_all.png"

rule plot:
    """그래프 그리기."""
    input:
        "../temp/wc_all.csv"
    output:
        "../temp/wc_all.png"
    script:
        "plot.py"
```

데이터 분석가 (u2) 는 데이터 엔지니어 (u1) 의 워크 플로우를 `include` 하고, 그것의 최종 타겟인 `../temp/wc_all.csv` 를 입력으로 `plot` 규칙을 정의하고 있다. 기본 규칙은 include 된 파일의 것은 무시하고, 현재 워크플로우의 첫 규칙인 `u2all` 이 선택한다.

이런 식으로, 중복 작업없이 다른 사람이 작업한 것을 가져다 사용할 수 있을 것이다.

## 기타 팁들

### 강제로 모든 규칙 실행하기

`--forceall` 줄여서 `-F` 옵션으로 유/무효 여부에 관계없이 모든 규칙을 실행할 수 있다.

```
$ snakemake -j -F
```

### 빌드 결과물 모두 지우기

예제처럼 하나의 디렉토리에 중간 결과물과 최종 타겟이 모두 저장되는 경우는 디렉토리 자체를 지우면 되나, 산출물 디렉토리를 몇 개로 구분해서 사용하는 경우에는 번거로울 수 있다. 이때는 아래 명령으로 빌드시 생성된 모든 파일을 제거할 수 있다.

```
$ snakemake -j --delete-all-output
```

## 마무리

Snakemake 에는 여기에 설명하지 않은 많은 기능들이 있다. 특히 [서브 워크플로우](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#sub-workflows) 같은 기능을 잘 활용하면, 다른 사람이 작성한 워크플로우의 산출물을 손쉽게 이용할 수 있어 협업이 용이해진다. 공식 문서 및 검색을 통해 Snakemake 의 강력한 기능을 좀 더 알아보고 활용하도록 하자.
