import React from 'react';
import { ScoreTier } from '../../types';
import { confTier, TIER_PILL_CLASS } from '../../lib/evidence';

export const TierPill: React.FC<{ tier: ScoreTier | string; className?: string }> = ({
  tier,
  className = '',
}) => {
  const t = confTier(tier);
  return (
    <span
      className={`inline-block text-[10px] font-bold px-2 py-0.5 rounded-full border ${TIER_PILL_CLASS[t.tone]} ${className}`}
    >
      {t.label}
    </span>
  );
};
