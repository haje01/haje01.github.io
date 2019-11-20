---
layout: template
title: RLlib Models, Preprocessors, and Action Distributions
tags: ray study
---

## RLlib의 컴포넌트들 간 데이터 흐름도

![](/assets/2019-11-20-11-33-04.png)

* 컴포넌트들
	* Environment - 환경
	* Preprocessor & Filter - 관측 전처리
	* Model - 뉴럴넷
	* ActionDistribution - 모델의 출력을 해석해 다음 동작 결정
* 녹색 컴포넌트들은 유저가 정의한 커스텀 구현으로 대체 가능

## 기본 동작
### 내장된 모델과 전처리기
* RLlib은 단순 휴리스틱으로 기본 모델을 선택
	* Vision network - 이미지 관측
	* FC - 나머지 전부
* `model`  컨피그 키로 설정
	* 관측이 커스텀 사이즈이면 `conv_filters`를 아래처럼 설정해야 할 것
	`"model": {"dim": 42, "conv_filters": [[16, [4, 4], 2], [32, [4, 4], 2], [512, [11, 11], 1]]}`
* `"model": {"use_lstm": true}`을 주면 모델 출력이 LSTM으로 처리될 것.
	* RLlib은 PG 알고리즘(A3C, PPO, PG, IMPALA)에 RNN을 지원
* RLlib은 환경의 관측 공간에 따라 내장 전처리기를 선택
	* 이산 관측은 OHE
	* Atari는 다운스케일 (DeepMind 전처리기 사용)
	* Tuple & Dict는 펼쳐짐

### TensorFlow 모델
* 커스텀 TF모델은 `TFModelV2`를 상속받아야 함.
	* 기존 `rllib.model.Model`을 대체
* `value_function`을 재정의해 커스텀 밸류 메소드를 구현
* `custom_loss`를 통해 지도(Supervised) / 자기-지도(Self-Supervised) 손실을 추가

#### 순환(Recurrent) 모델
* `use_lstm` 대신 커스텀 순환 모델 정의 가능
* `RecurrentTFModelV2`를 상속

### PyTorch 모델
* 같은 식으로 PyTorch 기반 알고리즘들(A2C, PG, QMIX)를 위한한 커스텀 PyTorch 모델을 추가 가능
* `TorchModelV2` 를 상속

## 커스텀 전처리기
* `ray.rllib.models.preprocessors.Preprocessor`을 상속받아 커스텀 전처리 구현

## 커스텀 동작 분포(Custom Action Distributions)
* 모델과 전처리기처럼 커스텀 동작 분포 클래스를 명시 가능
* `model`에 동작 분포 클래스의 레퍼런스를 건네어 사용
* **자동회귀 동작 출력** 구현에 주로 사용

```python
import ray
import ray.rllib.agents.ppo as ppo
from ray.rllib.models import ModelCatalog
from ray.rllib.models.preprocessors import Preprocessor

class MyActionDist(ActionDistribution):
    @staticmethod
    def required_model_output_shape(action_space, model_config):
        return 7  # 모델 출력 특성 벡터 크기를 조절

    def __init__(self, inputs, model):
        super(MyActionDist, self).__init__(inputs, model)
        assert model.num_outputs == 7

    def sample(self): ...
    def logp(self, actions): ...
    def entropy(self): ...

# 모델 카탈로그에 커스텀 동작 분포 등록
ModelCatalog.register_custom_action_dist("my_dist", MyActionDist)

ray.init()
trainer = ppo.PPOTrainer(env="CartPole-v0", config={
    "model": {
        # 커스텀 동작 분포 설정
        "custom_action_dist": "my_dist",
    },
})
```

## 지도 모델 손실(Supervised Model Losses)
* 모든 RLlib 알고리즘은 커스텀 모델을 통해 지도 손실 활용가능
	* 전문가 경험에 대한 모방학습 손실
	* 자기-지도 오토인코더 손실
* 이 손실은 정책 평가 입력이나 오프라인 저장소에서 읽은 데이터에 대해 정의 가능

### TensorFlow
* `custom_loss()`를 재정의해 커스텀 TF모델에 지도 손실을 추가
	* 여기에서 알고리즘의 기존 정책 손실을 받고 자신의 지도 손실을 추가하여 반환
* 오프라인 데이터에 대한 [Imitation Loss가 추가된 CartPole](https://github.com/ray-project/ray/blob/master/rllib/examples/custom_loss.py)예제

### PyTorch
* 정책 정의의 손실을 직접 변경 가능

## 가변 길이/모수화된 동작 공간(Variable-length / Parametric Action Spaces)
* 다음과 같은 환경에서 커스텀 모델 이용 가능
	* 스텝 당 유효 동작이 다양
	* 유효한 동작의 수가 아주 큼 -> 즉, 특정 관측에 따라 가능한 동작이 조건화
	* DQN이나 PG 계열 알고리즘에 사용가능
* 핵심 아이디어는 **동작의 의미는 관측에 완전히 조건화된다**는 것

1) 환경은 매 스텝마다 마스크와(또는) 유효한 동작의 임베딩을 관측의 일부로 함께 보내주어야 한다.
	* 배치가 가능하도록 동작의 수는 1에서 max 까지 변할 수 있다.

```python
class MyParamActionEnv(gym.Env):
    def __init__(self, max_avail_actions):
        self.action_space = Discrete(max_avail_actions)
        self.observation_space = Dict({
            "action_mask": Box(0, 1, shape=(max_avail_actions, )),
            "avail_actions": Box(-1, 1, shape=(max_avail_actions, action_embedding_sz)),
            "real_obs": ...,
        })
```

2) 관측의 `action_mask`와 `avail_actions` 부분을 **해석할 수 있는 커스텀 모델**을 정의. 모델은 네트웍의 출력과 각 동작 임베딩의 내적으로 액션 로짓을 계산. 무효한 동작은 확률을 0으로 스케일링하여 제외(Mask out)될 수 있음.
```python
class ParametricActionsModel(TFModelV2):
    def __init__(self,
                 obs_space,
                 action_space,
                 num_outputs,
                 model_config,
                 name,
                 true_obs_shape=(4,),
                 action_embed_size=2):
        super(ParametricActionsModel, self).__init__(
            obs_space, action_space, num_outputs, model_config, name)
        self.action_embed_model = FullyConnectedNetwork(...)

    def forward(self, input_dict, state, seq_lens):
        # 관측에서 유효한 동작 텐서를 추출.
        avail_actions = input_dict["obs"]["avail_actions"]
        action_mask = input_dict["obs"]["action_mask"]

        # 예측된 동작 임베딩을 계산
        action_embed, _ = self.action_embed_model({
            "obs": input_dict["obs"]["cart"]
        })

        # 임베딩을 [BATCH, 1, EMBED_SIZE]으로 확장. 유효한 동작 텐서는
        # [BATCH, MAX_ACTIONS, EMBED_SIZE] 형태인 것에 주목.
        intent_vector = tf.expand_dims(action_embed, 1)

        # 내적을 배치화 => 로짓의 형태는 [BATCH, MAX_ACTIONS].
        action_logits = tf.reduce_sum(avail_actions * intent_vector, axis=2)

        # 무효한 동작을 제외 (안정성을 위해 tf.float32.min 을 이용)
        # 가능한 동작 -> 0, 불가능한 동작 -> float32.min 으로 변경
        inf_mask = tf.maximum(tf.log(action_mask), tf.float32.min)

        # 불가능한 동작의 logit은 float32.min에 가깝게 된다.
        return action_logits + inf_mask, state
```

* 용도에 따라 마스킹 / 임베딩 중 하나만 하거나, 둘 다 할 수도 있다.
	* [parametric_action_cartpole.py](https://github.com/ray-project/ray/blob/master/rllib/examples/parametric_action_cartpole.py)예를 참고
	* 마스킹을 하면 모델 출력에  `tf.float32.min` 값이 나오기에, 모든 알고리즘과 동작 않을 수 있음.
		* 예로, 알고리즘이 `tf.float32.min` 값을 잘 처리하지 못하면 크래쉬가 발생
	* 카트폴 예는 DQN(`hiddens=[]` 설정 필요) 과 PPO(이동 평균을 끄고, `vf_share_layers=True`로 설정)및 몇가지 알고리즘과 잘 동작했다.
	* 모든 알고리즘이 모수화된 동작과 잘 동작하지는 않음

![](/assets/2019-11-20-11-33-18.png)

## 자동회귀 동작 분포(Autoregressive Action Distributions)
* 다중 요소가 있는 동작 공간(예: `Tuple(a1, a2)`)에서, `a2`가 표집된 `a1`의 값에 조건화되기를 바랄 수 있다.
	* 즉, `a2_sampled ~ P(a2 | a1_sampled, obs)`
	* 일반적으로 `a1`과 `a2`는 독립적으로 표집되는데, 정책의 표현력을 감소시킨다.
* 이를 위해, **자동회귀 패턴을 구현하는 커스텀 모델** 및 **이 모델을 활용하는 커스텀 동작 분포 클래스**가 필요하다.
* [autoregressive_action_dist.py](https://github.com/ray-project/ray/blob/master/rllib/examples/autoregressive_action_dist.py) 예는 단순 이진 동작 공간에서 구현을 보여줌.
* 더 복잡한 공간에서는, [MADE](https://arxiv.org/abs/1502.03509) 같은 더 효율적인 알고리즘이 필요.
* N-파트 동작은 모델의 N 전방 패스가 필요하나, **동작의 로그 확률 계산은 단일 패스에 가능**한 것에 주목

```python
class BinaryAutoregressiveOutput(ActionDistribution):
    """동작 분포 P(a1, a2) = P(a1) * P(a2 | a1)"""

    @staticmethod
    def required_model_output_shape(self, model_config):
        return 16  # 모델 출력 특성 벡터 크기를 조절

    def sample(self):
        # 먼저 a1을 표집
        a1_dist = self._a1_distribution()
        a1 = a1_dist.sample()

        # a1에 조건화된 a2를 표집
        a2_dist = self._a2_distribution(a1)
        a2 = a2_dist.sample()

        # 동작 튜플을 반환
        return TupleActions([a1, a2])

    def logp(self, actions):
        a1, a2 = actions[:, 0], actions[:, 1]
        a1_vec = tf.expand_dims(tf.cast(a1, tf.float32), 1)
        a1_logits, a2_logits = self.model.action_model([self.inputs, a1_vec])
        return (Categorical(a1_logits, None).logp(a1) + Categorical(
            a2_logits, None).logp(a2))

    def _a1_distribution(self):
        BATCH = tf.shape(self.inputs)[0]
        a1_logits, _ = self.model.action_model(
            [self.inputs, tf.zeros((BATCH, 1))])
        a1_dist = Categorical(a1_logits, None)
        return a1_dist

    def _a2_distribution(self, a1):
        a1_vec = tf.expand_dims(tf.cast(a1, tf.float32), 1)
        _, a2_logits = self.model.action_model([self.inputs, a1_vec])
        a2_dist = Categorical(a2_logits, None)
        return a2_dist
```

```python
class AutoregressiveActionsModel(TFModelV2):
    """위 코드에서 필요한 `.action_model` 브랜치를 구현."""

    def __init__(self, obs_space, action_space, num_outputs, model_config,
                 name):
        super(AutoregressiveActionsModel, self).__init__(
            obs_space, action_space, num_outputs, model_config, name)
        if action_space != Tuple([Discrete(2), Discrete(2)]):
            raise ValueError(
                "This model only supports the [2, 2] action space")

        # 입력
        obs_input = tf.keras.layers.Input(
            shape=obs_space.shape, name="obs_input")
        a1_input = tf.keras.layers.Input(shape=(1, ), name="a1_input")
        ctx_input = tf.keras.layers.Input(
            shape=(num_outputs, ), name="ctx_input")

        # 모델의 출력 (보통 'logits', 그러나 자동회귀 분포에서 이것은 관측을 인코딩하는 맥락/특성 레이어에 가까움)
        context = tf.keras.layers.Dense(
            num_outputs,
            name="hidden",
            activation=tf.nn.tanh,
            kernel_initializer=normc_initializer(1.0))(obs_input)

        # P(a1 | obs)
        a1_logits = tf.keras.layers.Dense(
            2,
            name="a1_logits",
            activation=None,
            kernel_initializer=normc_initializer(0.01))(ctx_input)

        # P(a2 | a1)
        # 주의: 보통은 다음처럼 P(a2 | a1, obs)을 구현하기를 원할 것:
        # a2_context = tf.keras.layers.Concatenate(axis=1)(
        #     [ctx_input, a1_input])
        a2_context = a1_input
        a2_hidden = tf.keras.layers.Dense(
            16,
            name="a2_hidden",
            activation=tf.nn.tanh,
            kernel_initializer=normc_initializer(1.0))(a2_context)
        a2_logits = tf.keras.layers.Dense(
            2,
            name="a2_logits",
            activation=None,
            kernel_initializer=normc_initializer(0.01))(a2_hidden)

        # 베이스 레이어
        self.base_model = tf.keras.Model(obs_input, context)
        self.register_variables(self.base_model.variables)
        self.base_model.summary()

        # 자동회귀 동작 샘플러
        self.action_model = tf.keras.Model([ctx_input, a1_input],
                                           [a1_logits, a2_logits])
        self.action_model.summary()
        self.register_variables(self.action_model.variables)
```

## 참고 링크
* https://ray.readthedocs.io/en/latest/rllib-models.html#autoregressive-action-distributions
