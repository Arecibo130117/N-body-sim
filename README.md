# N-body Simulator (Three.js + Vite)

이 프로젝트는 **Barnes-Hut 중력 근사 + 충돌/파편화 + 재응집(merge) + 발열/연소 VFX** 를 포함한
게임형 N-body 시뮬레이터 예제입니다.

> 주의: "행성 찌그러짐/발열/연소"는 **실시간 게임용 근사**(shader displacement + crater field + particle FX)입니다.
> 실제 유체/탄성체(FEM/MPM/SPH) 수준의 물리 정확도를 목표로 하진 않습니다.

## 실행

```bash
npm install
npm run dev
```

빌드:

```bash
npm run build
npm run preview
```

## 조작

- 좌클릭: 천체 선택
- 드래그: 카메라 회전
- 휠: 줌
- **SHIFT + 좌클릭**: 마우스 위치(원점 평면)에 천체 투하(질량/속도는 우측 패널)

## 구현 포인트

- `src/sim/Simulation.js` : 적분기(leapfrog), Barnes-Hut 옵션, merge/fragment/reaccrete
- `src/sim/collisions.js` : 상대속도/특정충돌에너지 기반으로 bounce/merge/fragment 결정
- `src/render/shaders.js` : crater field + damage wobble + thermal emissive(용암 균열 느낌)
- `src/render/Particles.js` : 파편/먼지/불꽃(가산 블렌딩) 파티클
